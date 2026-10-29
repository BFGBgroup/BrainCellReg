"""Database step 2 (draw-on-demand plots): the extra tables the plotting functions read, under
/public/home/yucong/scbrain/db/tables (+ db/demo for the sample pages), with checks appended to db/tables/CHECK2.txt.
  expression       long table gene x context x sample (deconvolved expression), sorted by gene for fast lookup
  expr_summary     gene x context: n, mean, sd, quartiles
  geno_box         box statistics of expression by genotype (0/1/2 copies of the effect allele) for every lead eSNP of an eGene
                   and every independent eSNP pair. Summary level only: no per-person genotype leaves the server.
  demo/            nominal cis results across the 110 contexts, GWAS P and 1000G-EUR LD for a few demo genes (locus plots)
Usage: python s26_db_plot_tables.py expr | geno | demo GENE1,GENE2,...
"""
import os, sys, glob
import numpy as np
import pandas as pd
import duckdb

B = '/public/home/yucong/scbrain'
IN = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data'
T = B + '/db/tables'
D = B + '/db/demo'
os.makedirs(D, exist_ok=True)
con = duckdb.connect()
con.execute("SET threads TO 8")
con.execute(f"SET temp_directory = '{B}/tmp/duck'")
con.execute("SET memory_limit = '48GB'")


def log_check(lines):
    with open(T + '/CHECK2.txt', 'a') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines), flush=True)


def read_ctx(f):
    ctx = os.path.basename(f).replace('.bed', '')
    d = pd.read_csv(f, sep='\t').drop(columns=['#chr', 'start', 'end']).set_index('gene-id')
    return ctx, d


def part_expr():
    files = sorted(glob.glob(IN + '/expression/110/*.bed'))
    tmp = B + '/tmp/expr_parts'
    os.makedirs(tmp, exist_ok=True)
    n_rows, n_files = 0, 0
    for f in files:
        ctx, d = read_ctx(f)
        x = d.reset_index().melt(id_vars='gene-id', var_name='sample', value_name='expr').rename(columns={'gene-id': 'gene'})
        x['sample'] = x['sample'].str.replace('-', '_', regex=False)
        x.insert(1, 'ctx', ctx)
        x['expr'] = x.expr.astype('float32')
        con.register('x', x)
        con.execute(f"COPY x TO '{tmp}/{ctx}.parquet' (FORMAT PARQUET)")
        con.unregister('x')
        n_rows += len(x); n_files += 1
    con.execute(f"""COPY (SELECT * FROM read_parquet('{tmp}/*.parquet') ORDER BY gene, ctx, sample)
                    TO '{T}/expression.parquet' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)""")
    con.execute(f"""COPY (SELECT gene, ctx, count(*) n, avg(expr) mean, stddev_samp(expr) sd,
                           quantile_cont(expr, 0.25) q1, median(expr) median, quantile_cont(expr, 0.75) q3
                          FROM read_parquet('{T}/expression.parquet') GROUP BY gene, ctx ORDER BY gene, ctx)
                    TO '{T}/expr_summary.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)""")
    e = con.execute(f"SELECT count(*), count(DISTINCT ctx), count(DISTINCT gene), count(DISTINCT sample) FROM '{T}/expression.parquet'").fetchone()
    s = con.execute(f"SELECT count(*) FROM '{T}/expr_summary.parquet'").fetchone()[0]
    # every (ctx, sample) of the expression matrices should be one of the 1,231 samples (region + individual)
    miss = con.execute(f"""SELECT count(*) FROM (SELECT DISTINCT split_part(ctx, '.', 1) region, sample FROM '{T}/expression.parquet') e
                           ANTI JOIN '{T}/samples.parquet' s ON e.region = s.region AND e.sample = s.individual""").fetchone()[0]
    log_check([f"{'OK ' if e[1] == 110 else 'BAD'} expression contexts: {e[1]} (110); files {n_files}",
               f"INFO expression rows {e[0]:,} (melted {n_rows:,}); genes {e[2]:,}; individuals {e[3]}",
               f"{'OK ' if miss == 0 else 'BAD'} expression (region, individual) not in samples table: {miss}",
               f"INFO expr_summary rows {s:,}"])


def read_bed_rows(idx, n_ind):
    bpr = (n_ind + 3) // 4
    mm = np.memmap(IN + '/genotype/out.bed', dtype=np.uint8, mode='r')
    assert mm[0] == 0x6c and mm[1] == 0x1b and mm[2] == 1
    raw = np.stack([mm[3 + i * bpr: 3 + (i + 1) * bpr] for i in idx])
    codes = np.stack([(raw >> (2 * k)) & 3 for k in range(4)], axis=2).reshape(len(idx), -1)[:, :n_ind]
    # plink 1 code: 0 = hom allele1, 1 = missing, 2 = het, 3 = hom allele2 -> copies of bim allele 2
    g = np.choose(codes, [0.0, np.nan, 1.0, 2.0])
    return g


def part_geno():
    fam = pd.read_csv(IN + '/genotype/out.fam', sep=r'\s+', header=None, usecols=[1])[1].str.replace('-', '_', regex=False).tolist()
    bim = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[1], names=['rsid'])
    bim['i'] = np.arange(len(bim))
    bim = bim.drop_duplicates('rsid').set_index('rsid').i
    es = con.execute(f"SELECT ctx, gene, lead_rsid rsid, slope, qval, effect_allele, other_allele, allele_status FROM '{T}/egene_stats.parquet' WHERE is_egene").df()
    es['lead'] = True
    pr = con.execute(f"SELECT ctx, gene, rsid, slope, effect_allele, other_allele, allele_status FROM '{T}/esnp_pairs.parquet' WHERE independent").df()
    pairs = pd.concat([es, pr], ignore_index=True).drop_duplicates(['ctx', 'gene', 'rsid'])
    pairs['lead'] = pairs.lead.eq(True)
    pairs = pairs[pairs.rsid.isin(bim.index)]
    snps = pairs.rsid.unique()
    G = read_bed_rows(bim.loc[snps].values, len(fam))
    gi = {s: k for k, s in enumerate(snps)}
    gdf_cols = {iid: j for j, iid in enumerate(fam)}
    out, chk = [], []
    for ctx, sub in pairs.groupby('ctx'):
        _, d = read_ctx(f'{IN}/expression/110/{ctx}.bed')
        d.columns = [c.replace('-', '_') for c in d.columns]
        cols = [gdf_cols[c] for c in d.columns]
        for r in sub.itertuples(index=False):
            if r.gene not in d.index:
                continue
            y = d.loc[r.gene].values.astype(float)
            g = G[gi[r.rsid], cols]
            ok = ~np.isnan(g)
            if ok.sum() < 10:
                continue
            gg, yy = g[ok], y[ok]
            if r.lead and r.qval < 1e-6 and gg.std() > 0:
                chk.append(np.sign(np.corrcoef(gg, yy)[0, 1]) == np.sign(r.slope))
            for k in (0, 1, 2):
                v = yy[gg == k]
                if len(v) == 0:
                    continue
                q1, med, q3 = np.percentile(v, [25, 50, 75])
                iqr = q3 - q1
                lo = v[v >= q1 - 1.5 * iqr].min(); hi = v[v <= q3 + 1.5 * iqr].max()
                out.append((ctx, r.gene, r.rsid, r.lead, k, len(v), float(v.mean()), float(q1), float(med), float(q3), float(lo), float(hi)))
    agree = float(np.mean(chk))
    flip = agree < 0.5
    gb = pd.DataFrame(out, columns=['ctx', 'gene', 'rsid', 'lead', 'g_bim_a2', 'n', 'mean', 'q1', 'median', 'q3', 'whisker_lo', 'whisker_hi'])
    # tensorQTL slopes are per copy of the counted allele (= effect_allele of the allele table); express groups the same way
    gb['n_effect'] = (2 - gb.g_bim_a2) if flip else gb.g_bim_a2
    gb = gb.drop(columns='g_bim_a2')
    gb = gb.merge(pairs[['ctx', 'gene', 'rsid', 'effect_allele', 'other_allele', 'allele_status']], on=['ctx', 'gene', 'rsid'], how='left')
    con.register('gb', gb)
    con.execute(f"COPY (SELECT * FROM gb ORDER BY gene, rsid, ctx, n_effect) TO '{T}/geno_box.parquet' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)")
    npair = gb[['ctx', 'gene', 'rsid']].drop_duplicates().shape[0]
    log_check([f"INFO geno_box pairs requested {len(pairs):,}, written {npair:,}; rows {len(gb):,}",
               f"{'OK ' if (agree > 0.95 or agree < 0.05) else 'BAD'} genotype coding vs tensorQTL slope sign (lead eSNPs q<1e-6, n={len(chk)}): "
               f"bim allele2 agrees {agree:.4f} -> counted allele = bim allele{'1' if flip else '2'}",
               f"INFO groups with n<3 (shown as points, not boxes): {int((gb.n < 3).sum()):,}"])


def part_demo(genes_csv):
    g = con.execute(f"SELECT gene, symbol, chr, start, \"end\", strand FROM '{T}/genes.parquet'").df()
    syms = genes_csv.split(',')
    g = g[g.symbol.isin(syms) | g.gene.isin(syms)]
    bimpos = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[1, 3], names=['rsid', 'pos38']).drop_duplicates('rsid')
    gw = {'PD': 'GCST009325', 'AD': 'GCST013196', 'SCZ': 'GCST90128471', 'MDD': 'GCST90726344'}
    bf = B + '/local_handoff/scbrain1/g1000_eur/g1000_eur'
    for r in g.itertuples(index=False):
        files = sorted(glob.glob(f'{IN}/run/new/*.cis_qtl_pairs.chr{r.chr}.parquet'))
        parts = []
        for f in files:
            ctx = os.path.basename(f).split('.cis_qtl_pairs')[0]
            x = con.execute(f"SELECT variant_id rsid, pval_nominal p, slope, slope_se FROM read_parquet('{f}') WHERE phenotype_id = '{r.gene}'").df()
            x.insert(0, 'ctx', ctx)
            parts.append(x)
        nom = pd.concat(parts, ignore_index=True).merge(bimpos, on='rsid', how='left')
        con.register('nom', nom)
        con.execute(f"COPY nom TO '{D}/{r.symbol}_nominal.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)")
        con.unregister('nom')
        snps = nom.rsid.drop_duplicates()
        gws = []
        for dz, acc in gw.items():
            x = pd.read_csv(f'{B}/gwas10/{acc}-smr.tsv', sep='\t', usecols=['SNP', 'p'])
            x = x[x.SNP.isin(snps)].rename(columns={'SNP': 'rsid'})
            x.insert(0, 'disease', dz)
            gws.append(x)
        gwdf = pd.concat(gws, ignore_index=True).merge(bimpos, on='rsid', how='left')
        con.register('gw', gwdf)
        con.execute(f"COPY gw TO '{D}/{r.symbol}_gwas.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)")
        con.unregister('gw')
        # LD (1000G EUR, r2) with every context's lead variant and each GWAS's top variant in the window
        idx = set(nom.loc[nom.groupby('ctx').p.idxmin(), 'rsid'])
        for dz, x in gwdf.groupby('disease'):
            idx.add(x.loc[x.p.idxmin(), 'rsid'])
        lst = f'{D}/{r.symbol}_ldidx.txt'
        snpl = f'{D}/{r.symbol}_snps.txt'
        open(lst, 'w').write('\n'.join(sorted(idx)) + '\n')
        snps.to_csv(snpl, index=False, header=False)
        os.system(f"{B}/bin/plink --bfile {bf} --extract {snpl} --r2 --ld-snp-list {lst} --ld-window-kb 2000 --ld-window 999999 "
                  f"--ld-window-r2 0 --out {D}/{r.symbol}_ld > /dev/null 2>&1")
        ld = pd.read_csv(f'{D}/{r.symbol}_ld.ld', sep=r'\s+', usecols=['SNP_A', 'SNP_B', 'R2'])
        ld = ld.rename(columns={'SNP_A': 'index_rsid', 'SNP_B': 'rsid', 'R2': 'r2'})
        con.register('ld', ld)
        con.execute(f"COPY ld TO '{D}/{r.symbol}_ld.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)")
        con.unregister('ld')
        log_check([f"INFO demo {r.symbol} ({r.gene}, chr{r.chr}): nominal rows {len(nom):,} in {nom.ctx.nunique()} contexts, "
                   f"variants {len(snps):,}; GWAS rows {len(gwdf):,}; LD index {len(idx)}, LD rows {len(ld):,}"])


if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'expr':
        part_expr()
    elif what == 'geno':
        part_geno()
    elif what == 'demo':
        part_demo(sys.argv[2])
    print(what.upper() + '_DONE', flush=True)
