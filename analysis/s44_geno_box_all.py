"""Box statistics by genotype for EVERY significant eSNP-gene pair (PI #118, 10-17: rs156429 had boxes in 6 of its 34 significant
contexts, because s26 kept only lead and independent pairs). Same statistics and allele coding as s26 part_geno, vectorised per context.
Writes /public/home/yucong/scbrain/db/geno_box_all.parquet (outside tables/, so the site's current table is untouched) and
db/CHECK3.txt: pair counts, allele coding check, and an exact comparison with the 214,551 pairs of the old geno_box.
Summary level only: no per-person genotype leaves the server.
"""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s26_db_plot_tables import read_ctx, IN, T, B, con

OUT = B + '/db/geno_box_all.parquet'
CHK = B + '/db/CHECK3.txt'
CHUNK = 20000


def log(lines):
    with open(CHK, 'a') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines), flush=True)


def bed_rows(idx, n_ind, n_snp):
    """int8 copies of bim allele 2 (-1 = missing) for the bim rows idx."""
    bpr = (n_ind + 3) // 4
    mm = np.memmap(IN + '/genotype/out.bed', dtype=np.uint8, mode='r')
    assert mm[0] == 0x6c and mm[1] == 0x1b and mm[2] == 1
    body = mm[3:3 + n_snp * bpr].reshape(n_snp, bpr)
    lut = np.array([0, -1, 1, 2], dtype=np.int8)  # plink 1: 0 hom a1, 1 missing, 2 het, 3 hom a2
    out = np.empty((len(idx), n_ind), dtype=np.int8)
    for s in range(0, len(idx), 200000):
        raw = np.asarray(body[idx[s:s + 200000]])
        codes = np.stack([(raw >> (2 * k)) & 3 for k in range(4)], axis=2).reshape(len(raw), -1)[:, :n_ind]
        out[s:s + 200000] = lut[codes]
    return out


def box_stats(Y, Gm):
    """Per row and genotype k: n, mean, q1, median, q3, whiskers (1.5 IQR rule, as s26)."""
    res = {}
    for k in (0, 1, 2):
        V = np.where(Gm == k, Y, np.nan)
        n = np.sum(Gm == k, axis=1)
        with np.errstate(all='ignore'):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                q1, med, q3 = np.nanpercentile(V, [25, 50, 75], axis=1)
                mean = np.nanmean(V, axis=1)
                iqr = q3 - q1
                lo = np.nanmin(np.where(V >= (q1 - 1.5 * iqr)[:, None], V, np.nan), axis=1)
                hi = np.nanmax(np.where(V <= (q3 + 1.5 * iqr)[:, None], V, np.nan), axis=1)
        res[k] = (n, mean, q1, med, q3, lo, hi)
    return res


def main():
    open(CHK, 'w').close()
    fam = pd.read_csv(IN + '/genotype/out.fam', sep=r'\s+', header=None, usecols=[1])[1].str.replace('-', '_', regex=False).tolist()
    bim = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[1], names=['rsid'])
    n_snp = len(bim)
    bim['i'] = np.arange(n_snp)
    bim = bim.drop_duplicates('rsid').set_index('rsid').i
    es = con.execute(f"SELECT ctx, gene, lead_rsid rsid, slope, qval, effect_allele, other_allele, allele_status FROM '{T}/egene_stats.parquet' WHERE is_egene").df()
    es['lead'] = True
    pr = con.execute(f"SELECT ctx, gene, rsid, slope, effect_allele, other_allele, allele_status FROM '{T}/esnp_pairs.parquet'").df()
    n_sig = len(pr)
    pairs = pd.concat([es, pr], ignore_index=True).drop_duplicates(['ctx', 'gene', 'rsid'])
    pairs['lead'] = pairs.lead.eq(True)
    n_req = len(pairs)
    pairs = pairs[pairs.rsid.isin(bim.index)].reset_index(drop=True)
    snps = pairs.rsid.unique()
    G = bed_rows(bim.loc[snps].values, len(fam), n_snp)
    gi = pd.Series(np.arange(len(snps)), index=snps)
    fcol = {iid: j for j, iid in enumerate(fam)}
    log([f"INFO esnp_pairs {n_sig:,}; with lead eSNPs {n_req:,}; in genotype bim {len(pairs):,}; variants {len(snps):,}; individuals {len(fam)}"])
    frames, sgn = [], []
    for ctx, sub in pairs.groupby('ctx'):
        _, d = read_ctx(f'{IN}/expression/110/{ctx}.bed')
        d.columns = [c.replace('-', '_') for c in d.columns]
        cols = np.array([fcol[c] for c in d.columns])
        sub = sub[sub.gene.isin(d.index)]
        for s in range(0, len(sub), CHUNK):
            c = sub.iloc[s:s + CHUNK]
            Y = d.loc[c.gene.values].values.astype(float)
            Gm = G[gi.loc[c.rsid.values].values][:, cols]
            ok = (Gm >= 0).sum(axis=1) >= 10
            c, Y, Gm = c[ok], Y[ok], Gm[ok]
            # allele coding check on strong lead pairs (as s26)
            m = (c.lead & (c.qval < 1e-6)).values
            if m.any():
                g = np.where(Gm[m] >= 0, Gm[m], np.nan).astype(float); y = np.where(Gm[m] >= 0, Y[m], np.nan)
                cov = np.nanmean((g - np.nanmean(g, 1, keepdims=True)) * (y - np.nanmean(y, 1, keepdims=True)), 1)
                sl = c.slope.values[m]
                keep = np.nanstd(g, 1) > 0
                sgn.extend((np.sign(cov[keep]) == np.sign(sl[keep])).tolist())
            st = box_stats(Y, Gm)
            for k in (0, 1, 2):
                n, mean, q1, med, q3, lo, hi = st[k]
                h = n > 0
                frames.append(pd.DataFrame({'ctx': ctx, 'gene': c.gene.values[h], 'rsid': c.rsid.values[h], 'lead': c.lead.values[h], 'g_bim_a2': k,
                                            'n': n[h], 'mean': mean[h], 'q1': q1[h], 'median': med[h], 'q3': q3[h], 'whisker_lo': lo[h], 'whisker_hi': hi[h]}))
        print(ctx, len(sub), flush=True)
    agree = float(np.mean(sgn))
    flip = agree < 0.5
    gb = pd.concat(frames, ignore_index=True)
    gb['n_effect'] = (2 - gb.g_bim_a2) if flip else gb.g_bim_a2
    gb = gb.drop(columns='g_bim_a2')
    gb = gb.merge(pairs[['ctx', 'gene', 'rsid', 'effect_allele', 'other_allele', 'allele_status']], on=['ctx', 'gene', 'rsid'], how='left')
    con.register('gb', gb)
    con.execute(f"COPY (SELECT * FROM gb ORDER BY gene, rsid, ctx, n_effect) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 50000)")
    npair = con.execute(f"SELECT count(*) FROM (SELECT DISTINCT ctx, gene, rsid FROM '{OUT}')").fetchone()[0]
    # exact comparison with the old table on the pairs both have
    cmp = con.execute(f"""SELECT count(*) n, max(abs(a.mean - b.mean)) dmean, max(abs(a.median - b.median)) dmed,
                                 max(abs(a.whisker_lo - b.whisker_lo)) dlo, max(abs(a.whisker_hi - b.whisker_hi)) dhi, sum((a.n <> b.n)::int) dn
                          FROM '{T}/geno_box.parquet' a JOIN '{OUT}' b USING (ctx, gene, rsid, n_effect)""").fetchone()
    old_rows = con.execute(f"SELECT count(*) FROM '{T}/geno_box.parquet'").fetchone()[0]
    rs = con.execute(f"SELECT count(DISTINCT ctx) FROM '{OUT}' WHERE rsid = 'rs156429' AND gene = (SELECT gene FROM '{T}/genes.parquet' WHERE symbol = 'GPNMB')").fetchone()[0]
    log([f"INFO geno_box_all pairs written {npair:,}; rows {len(gb):,}; file {os.path.getsize(OUT) / 1e6:.0f} MB",
         f"{'OK ' if (agree > 0.95 or agree < 0.05) else 'BAD'} allele coding vs slope sign (lead q<1e-6, n={len(sgn)}): bim allele2 agrees {agree:.4f} -> counted allele = bim allele{'1' if flip else '2'}",
         f"{'OK ' if cmp[0] == old_rows and cmp[5] == 0 and max(cmp[1:5]) < 1e-9 else 'BAD'} old geno_box rows matched {cmp[0]:,} of {old_rows:,}; n differ {cmp[5]}; "
         f"max |diff| mean {cmp[1]:.2e} median {cmp[2]:.2e} whiskers {cmp[3]:.2e}/{cmp[4]:.2e}",
         f"INFO rs156429-GPNMB contexts with boxes: {rs} (significant in 34)"])


if __name__ == '__main__':
    main()
    print('GENO_ALL_DONE', flush=True)
