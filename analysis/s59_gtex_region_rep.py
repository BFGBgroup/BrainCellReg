"""分析 A（规格 progress/分析规格_GTEx脑区复现与皮层下疾病基因_v1.md）：GTEx v8 13 个脑组织的 eGene lead 对，在本稿 110 个组合里的
复现（π1，λ=0.5）与方向一致率 → 10 脑区 × 13 组织矩阵、5 个对应对的按列 / 按行判定、置换检验。
  in   external/gtex_v8/QTD000*.permuted.tsv.gz, dataset_metadata_r7.tsv; nominal_web/<chr>/<ctx>.parquet; db/tables_allele_v3/variants.parquet
  out  gtex_rep/pairs.parquet（逐对查到的值）, gtex_rep/pi1_ctx.tsv, gtex_rep/matrix_pi1.tsv, matrix_concord.tsv, verdict.txt"""
import glob, os, sys
from multiprocessing import Pool
import duckdb
import numpy as np
import pandas as pd

B = '/public/home/yucong/scbrain'
G = B + '/external/gtex_v8'
OUT = B + '/gtex_rep'
os.makedirs(OUT, exist_ok=True)
LAMBDA = 0.5
REG = ['FCTX', 'TCTX', 'OCTX', 'HIPP', 'PUTM', 'THAL', 'SNIG', 'MEDU', 'WHMT', 'CRBL']
MATCH = {'FCTX': 'brain_frontal_cortex', 'HIPP': 'brain_hippocampus', 'PUTM': 'brain_putamen',
         'SNIG': 'brain_substantia_nigra', 'CRBL': 'brain_cerebellum'}
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def pi1(p):
    p = np.asarray(p, float)
    return float(np.clip(1 - (p > LAMBDA).mean() / (1 - LAMBDA), 0, 1)) if len(p) >= 20 else np.nan


meta = pd.read_csv(G + '/dataset_metadata_r7.tsv', sep='\t')
meta = meta[(meta.study_label == 'GTEx') & (meta.quant_method == 'ge') & meta.sample_group.str.startswith('brain_')]
leads = []
for r in meta.itertuples():
    d = pd.read_csv(f'{G}/{r.dataset_id}.permuted.tsv.gz', sep='\t')
    d = d[d.chromosome.astype(str).isin([str(i) for i in range(1, 23)])].copy()
    d['q'] = bh(d.p_beta.values)
    d = d[d.q < 0.05]
    v = d.variant.str.split('_', expand=True)
    d['ref'], d['alt'] = v[2], v[3]
    d['tissue'], d['n_gtex'] = r.sample_group, r.sample_size
    leads.append(d[['tissue', 'n_gtex', 'molecular_trait_id', 'chromosome', 'position', 'ref', 'alt', 'beta']])
leads = pd.concat(leads, ignore_index=True).rename(columns={'molecular_trait_id': 'gene', 'chromosome': 'chr', 'position': 'pos', 'beta': 'beta_g'})
leads['gene'] = leads.gene.str.replace(r'\..*$', '', regex=True)
print('GTEx eGenes (BH<0.05) per tissue'); print(leads.groupby('tissue').size().to_string(), flush=True)
var = duckdb.sql(f"SELECT rsid, a1, a2, effect_allele FROM '{B}/db/tables_allele_v3/variants.parquet'").df()


def one_chr(c):
    keys = leads[leads.chr.astype(str) == str(c)][['gene', 'pos']].drop_duplicates()
    if keys.empty:
        return None
    con = duckdb.connect(); con.register('k', keys)
    out = []
    for f in sorted(glob.glob(f'{B}/nominal_web/chr{c}/*.parquet')):
        ctx = os.path.basename(f)[:-8]
        x = con.execute(f"SELECT n.pos, n.phenotype_id AS gene, n.variant_id AS rsid, n.slope, n.mlog10p FROM read_parquet('{f}') n "
                        "JOIN k ON n.pos = k.pos AND n.phenotype_id = k.gene").df()
        x['ctx'] = ctx
        out.append(x)
    return pd.concat(out, ignore_index=True) if out else None


if __name__ == '__main__':
    with Pool(8) as pool:
        parts = [p for p in pool.map(one_chr, range(1, 23)) if p is not None]
    hit = pd.concat(parts, ignore_index=True)
    hit = hit.merge(var, on='rsid', how='left')
    m = leads.merge(hit, on=['gene', 'pos'], how='inner')
    # allele letters must agree (same or complementary strand); multi-allelic sites otherwise pair with the wrong variant
    same = ((m.a1 == m.ref) & (m.a2 == m.alt)) | ((m.a1 == m.alt) & (m.a2 == m.ref))
    comp = ((m.a1.map(COMP) == m.ref) & (m.a2.map(COMP) == m.alt)) | ((m.a1.map(COMP) == m.alt) & (m.a2.map(COMP) == m.ref))
    m = m[same | comp].copy()
    m['p'] = 10 ** (-m.mlog10p)
    ea = m.effect_allele
    m['sign_ours_alt'] = np.where(ea == m.alt, np.sign(m.slope), np.where(ea == m.ref, -np.sign(m.slope),
                         np.where(ea.map(COMP) == m.alt, np.sign(m.slope), np.where(ea.map(COMP) == m.ref, -np.sign(m.slope), np.nan))))
    pal = m.ref.map(COMP) == m.alt
    m.loc[pal & (m.effect_allele.isna()), 'sign_ours_alt'] = np.nan
    m['region'] = m.ctx.str.split('.').str[0]
    m.to_parquet(OUT + '/pairs.parquet', index=False)
    print('looked-up pairs', len(m), 'of GTEx leads', len(leads), flush=True)

    rows = []
    for (t, ctx), d in m.groupby(['tissue', 'ctx']):
        s = d[(d.p < 0.05) & d.sign_ours_alt.notna()]
        rows.append(dict(tissue=t, ctx=ctx, region=ctx.split('.')[0], n=len(d), pi1=pi1(d.p),
                         n_sig=len(s), concord=float((s.sign_ours_alt == np.sign(s.beta_g)).mean()) if len(s) >= 10 else np.nan))
    pc = pd.DataFrame(rows)
    pc.to_csv(OUT + '/pi1_ctx.tsv', sep='\t', index=False)
    M = pc.groupby(['region', 'tissue']).pi1.median().unstack().reindex(REG)
    C = pc.groupby(['region', 'tissue']).concord.median().unstack().reindex(REG)
    M.to_csv(OUT + '/matrix_pi1.tsv', sep='\t'); C.to_csv(OUT + '/matrix_concord.tsv', sep='\t')

    def hits(mat, byrow=False):
        h = {}
        for r, t in MATCH.items():
            h[r] = (mat.loc[r].idxmax() == t) if byrow else (mat[t].idxmax() == r)
        return h

    col = hits(M); row = hits(M, byrow=True)
    obs = sum(col.values())
    # shuffling the row labels: the label at position k gets the values of original row perm[k], so region r is the column max
    # of tissue t exactly when perm[index(r)] is the original argmax row of t
    rng = np.random.default_rng(1)
    amax = {t: REG.index(M[t].idxmax()) for t in MATCH.values()}
    perms = np.array([rng.permutation(len(REG)) for _ in range(100000)])
    stat = sum((perms[:, REG.index(r)] == amax[t]).astype(int) for r, t in MATCH.items())
    pperm = (int((stat >= obs).sum()) + 1) / 100001
    lines = ['π1 matrix (median over cell types; rows = this study regions, columns = GTEx tissues)', M.round(3).to_string(), '',
             'direction concordance matrix', C.round(3).to_string(), '',
             'by column (primary): ' + ', '.join(f'{r}:{"yes" if v else "no"}' for r, v in col.items()) + f'  -> {obs}/5, permutation P = {pperm:.4f}',
             'by row (secondary): ' + ', '.join(f'{r}:{"yes" if v else "no"}' for r, v in row.items()) + f'  -> {sum(row.values())}/5',
             'verdict: ' + ('PASS (>=4/5 and P<0.05)' if obs >= 4 and pperm < 0.05 else ('NOT IN PAPER (<3/5)' if obs < 3 else 'INCONCLUSIVE (3/5, or P>=0.05)'))]
    open(OUT + '/verdict.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print('S59_DONE')
