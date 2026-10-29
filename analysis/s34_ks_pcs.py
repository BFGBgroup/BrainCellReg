"""Fig 4d (PI #84 「enformer有20个PC，他这里enformer L2是什么？favor也有很多aPC」): KS of each of the 20 Enformer PCs,
eSNPs vs non-eSNPs, for WB, the 11 cell types and the 10 regions. Same sets, matching and direction rule as s19_ks_ctx.py
(L2 there = sqrt of the sum of the 20 squared PCs). Output: esnp_eg/results/ks_pcs_by_group.tsv"""
import os, glob, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy import stats

S = '/public/home/yucong/scbrain'
O = f'{S}/esnp_eg'
R = f'{O}/results'
REGIONS = ['CRBL', 'FCTX', 'HIPP', 'MEDU', 'OCTX', 'PUTM', 'SNIG', 'TCTX', 'THAL', 'WHMT']
CELLS = ['Astro', 'Endo', 'ExN', 'Gran', 'InN', 'Micro', 'NPC', 'Olig', 'OPC', 'Perc', 'Purk']
t0 = time.time()
var = pd.read_csv(f'{O}/sets/variants.tsv.gz', sep='\t', dtype={'chr': str})
m = np.load(f'{O}/sets/membership.npz', allow_pickle=True)
ctx = [str(x) for x in m['ctx']]
e3, n4, pre = m['e3'], m['n4'], m['pre']
eg = np.load(f'{O}/sets/membership_eg.npz')['eg']
nU = len(var)
G = {'WB': list(range(len(ctx)))}
for c in CELLS:
    G['cell:' + c] = [j for j, x in enumerate(ctx) if x.split('.')[1] == c]
for r in REGIONS:
    G['region:' + r] = [j for j, x in enumerate(ctx) if x.split('.')[0] == r]
MASK = {}
for g, J in G.items():
    E = e3[:, J].any(1)
    MASK[g] = ((e3[:, J] & eg[:, J]).any(1), n4[:, J].any(1) & ~pre[:, J].any(1) & ~E)
assert int(MASK['WB'][0].sum()) == 37412
var['idx'] = np.arange(nU)
PCS = [f'PC{i}' for i in range(20)]
enf = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['chrom', 'pos19', 'rsid', 'ref', 'alt'] + PCS, dtype={'ref': str, 'alt': str})
                 for f in sorted(glob.glob(f'{O}/annot/enformer/chr*.tsv')) if os.path.getsize(f) > 0], ignore_index=True)
enf = enf.merge(var[['idx', 'rsid', 'a1', 'a2']], on='rsid', how='inner')
enf['match'] = ((enf.ref == enf.a1) & (enf.alt == enf.a2)) | ((enf.ref == enf.a2) & (enf.alt == enf.a1))
enf = enf.sort_values(['idx', 'match'], ascending=[True, False]).drop_duplicates('idx')
FEAT = {}
for c in PCS:
    a = np.full(nU, np.nan, dtype=np.float32); a[enf['idx'].values] = enf[c].values.astype(np.float32); FEAT['enformer_' + c] = a
print(f'[{time.time() - t0:.0f}s] enformer variants', len(enf), flush=True)
del enf


def task(g):
    E, N = MASK[g]
    out = []
    for f, a in FEAT.items():
        x = a[E & ~np.isnan(a)]; y = a[N & ~np.isnan(a)]
        r = stats.ks_2samp(x, y, method='asymp')
        xs, ys = np.sort(x), np.sort(y)
        grid = np.concatenate([xs, ys])
        gap = np.searchsorted(ys, grid, 'right') / len(ys) - np.searchsorted(xs, grid, 'right') / len(xs)
        k = int(np.argmax(np.abs(gap)))
        out.append(dict(group=g, feature=f, n_esnp=len(x), n_non=len(y), median_esnp=float(np.median(x)), median_non=float(np.median(y)),
                        D=float(r.statistic), direction=('higher' if gap[k] > 0 else 'lower'), p=float(r.pvalue)))
    return out


with Pool(8) as pool:
    res = pool.map(task, list(G))
KS = pd.DataFrame([r for a in res for r in a])
KS.to_csv(f'{R}/ks_pcs_by_group.tsv', sep='\t', index=False)
print(KS[KS.group == 'WB'][['feature', 'D', 'direction', 'p']].to_string(), flush=True)
print('rows', len(KS), '| P<0.05', int((KS.p < 0.05).sum()))
print('KS_PCS_DONE')
