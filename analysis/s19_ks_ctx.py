"""Fig 4 v7 inputs (PI 10-14 #33): every brain region x cell type context, continuous scores as cumulative distributions + KS test.
Same sets as s4_stats.py with eset 'eg' (independent eSNPs of q<0.05 eGenes) vs non-eSNPs of the same group.
Writes to $SCB_OUT/results:
  ks_by_group.tsv     group (WB / region / cell / ctx), feature, n_esnp, n_non, median_esnp, median_non, D, direction, p, q (BH within feature x level)
  ecdf_by_group.tsv   group, feature, set (eSNP / non-eSNP), x, F   (51-point grid per feature, grid = pooled whole-brain 0.5-99.5 percentiles)
  tss_by_ctx.tsv      ctx, n_pairs, pct_10kb, pct_100kb, median_abs_dist  (eSNP-gene pairs of q<0.05 eGenes, before clumping)
"""
import os, sys, glob, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy import stats

S = '/public/home/yucong/scbrain'
O = os.environ.get('SCB_OUT', f'{S}/esnp_eg')
R = f'{O}/results'
NPROC = int(os.environ.get('NPROC', '8'))
REGIONS = ['CRBL', 'FCTX', 'HIPP', 'MEDU', 'OCTX', 'PUTM', 'SNIG', 'TCTX', 'THAL', 'WHMT']
CELLS = ['Astro', 'Endo', 'ExN', 'Gran', 'InN', 'Micro', 'NPC', 'Olig', 'OPC', 'Perc', 'Purk']
t0 = time.time()


def log(*a):
    print(f'[{time.time() - t0:7.1f}s]', *a, flush=True)


def bh(p):
    p = np.asarray(p, float); q = np.full(p.shape, np.nan); ok = ~np.isnan(p); n = ok.sum()
    if n == 0:
        return q
    o = np.argsort(p[ok]); r = p[ok][o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(r, 1); q[ok] = out
    return q


var = pd.read_csv(f'{O}/sets/variants.tsv.gz', sep='\t', dtype={'chr': str})
m = np.load(f'{O}/sets/membership.npz', allow_pickle=True)
ctx = [str(x) for x in m['ctx']]
e3, n4, pre = m['e3'], m['n4'], m['pre']
eg = np.load(f'{O}/sets/membership_eg.npz')['eg']
nU, nC = e3.shape
assert nU == len(var) and eg.shape == e3.shape
GROUPS = {'WB': list(range(nC))}
for c in CELLS:
    GROUPS['cell:' + c] = [j for j, x in enumerate(ctx) if x.split('.')[1] == c]
for r in REGIONS:
    GROUPS['region:' + r] = [j for j, x in enumerate(ctx) if x.split('.')[0] == r]
for j, x in enumerate(ctx):
    GROUPS['ctx:' + x] = [j]
ESET = os.environ.get('ESET', 'eg')            # eg = official set; lnc / pc = official set split by the biotype of the eGene
SUF = '' if ESET == 'eg' else '_' + ESET
if ESET != 'eg':
    eg = eg & m[ESET]
MASK = {}
for g, J in GROUPS.items():
    E = e3[:, J].any(1)
    Eg = (e3[:, J] & eg[:, J]).any(1)
    N = n4[:, J].any(1) & ~pre[:, J].any(1) & ~E
    MASK[g] = (Eg, N)
wb_e = int(MASK['WB'][0].sum())
log('eSNP set', ESET, '| universe', nU, 'contexts', nC, '| WB eSNPs', wb_e, '| WB non', int(MASK['WB'][1].sum()))
if ESET == 'eg':
    assert wb_e == 37412, 'whole-brain eSNP count differs from the official 37,412'
var['idx'] = np.arange(nU)

APC = ['apc_conservation', 'apc_conservation_v2', 'apc_epigenetics', 'apc_epigenetics_active',
       'apc_epigenetics_repressed', 'apc_epigenetics_transcription', 'apc_local_nucleotide_diversity',
       'apc_local_nucleotide_diversity_v2', 'apc_local_nucleotide_diversity_v3', 'apc_mappability',
       'apc_micro_rna', 'apc_mutation_density', 'apc_protein_function', 'apc_proximity_to_coding',
       'apc_proximity_to_coding_v2', 'apc_proximity_to_tsstes', 'apc_transcription_factor']
FCOLS = ['chr', 'pos', 'ref', 'alt'] + APC + ['rsid_f', 'category', 'linsight', 'cadd_phred']
USE = [a for a in APC if not a.endswith(('_v2', '_v3'))] + ['linsight', 'cadd_phred']
FEAT = {}
fl = sorted(glob.glob(f'{O}/annot/favor/chr*.tsv'))
fav = pd.concat([pd.read_csv(f, sep='\t', header=None, names=FCOLS, usecols=['chr', 'pos', 'ref', 'alt'] + USE, na_values=['\\N'],
                             dtype={'chr': str, 'ref': str, 'alt': str}) for f in fl if os.path.getsize(f) > 0], ignore_index=True)
fav = fav.merge(var[['idx', 'chr', 'pos', 'a1', 'a2']], on=['chr', 'pos'], how='inner')
fav['match'] = ((fav.ref == fav.a1) & (fav.alt == fav.a2)) | ((fav.ref == fav.a2) & (fav.alt == fav.a1))
fav = fav.sort_values(['idx', 'match'], ascending=[True, False]).drop_duplicates('idx')
log('favor matched variants', len(fav))
for c in USE:
    a = np.full(nU, np.nan, dtype=np.float32); a[fav['idx'].values] = fav[c].values.astype(np.float32); FEAT[c] = a
del fav
el = sorted(glob.glob(f'{O}/annot/enformer/chr*.tsv'))
PCS = [f'PC{i}' for i in range(20)]
enf = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['chrom', 'pos19', 'rsid', 'ref', 'alt'] + PCS, dtype={'ref': str, 'alt': str})
                 for f in el if os.path.getsize(f) > 0], ignore_index=True)
enf = enf.merge(var[['idx', 'rsid', 'a1', 'a2']], on='rsid', how='inner')
enf['match'] = ((enf.ref == enf.a1) & (enf.alt == enf.a2)) | ((enf.ref == enf.a2) & (enf.alt == enf.a1))
enf = enf.sort_values(['idx', 'match'], ascending=[True, False]).drop_duplicates('idx')
l2 = np.full(nU, np.nan, dtype=np.float32); l2[enf['idx'].values] = np.sqrt((enf[PCS].values ** 2).sum(1)).astype(np.float32)
FEAT['enformer_L2'] = l2
del enf
log('features', list(FEAT))

# grid per feature from the pooled whole-brain sets
E0, N0 = MASK['WB']
GRID = {}
for f, a in FEAT.items():
    v = a[(E0 | N0) & ~np.isnan(a)]
    lo, hi = np.percentile(v, [0.5, 99.5])
    GRID[f] = np.unique(np.linspace(lo, hi, 51))


def task(g):
    E, N = MASK[g]
    ks, ec = [], []
    for f, a in FEAT.items():
        x = a[E & ~np.isnan(a)]; y = a[N & ~np.isnan(a)]
        if len(x) < 5 or len(y) < 5:
            ks.append(dict(group=g, feature=f, n_esnp=len(x), n_non=len(y)))
            continue
        r = stats.ks_2samp(x, y, method='asymp')
        mx, my = float(np.median(x)), float(np.median(y))
        # direction: sign of the largest ECDF gap (non-eSNP CDF above eSNP CDF = eSNPs shifted to higher values)
        xs, ys = np.sort(x), np.sort(y)
        grid = np.concatenate([xs, ys])
        gap = np.searchsorted(ys, grid, 'right') / len(ys) - np.searchsorted(xs, grid, 'right') / len(xs)
        k = int(np.argmax(np.abs(gap)))
        ks.append(dict(group=g, feature=f, n_esnp=len(x), n_non=len(y), median_esnp=mx, median_non=my,
                       D=float(r.statistic), direction=('higher' if gap[k] > 0 else 'lower'), p=float(r.pvalue)))
        gx = GRID[f]
        ec.append(pd.DataFrame(dict(group=g, feature=f, set='eSNP', x=gx, F=np.searchsorted(xs, gx, 'right') / len(xs))))
        ec.append(pd.DataFrame(dict(group=g, feature=f, set='non-eSNP', x=gx, F=np.searchsorted(ys, gx, 'right') / len(ys))))
    return ks, ec


with Pool(NPROC) as pool:
    res = pool.map(task, list(GROUPS), chunksize=2)
KS = pd.DataFrame([r for a, _ in res for r in a])
KS['level'] = KS.group.str.split(':').str[0]
KS['q'] = np.nan
for (lev, f), i in KS.groupby(['level', 'feature']).groups.items():
    KS.loc[i, 'q'] = bh(KS.loc[i, 'p'].values)
KS.to_csv(f'{R}/ks_by_group{SUF}.tsv', sep='\t', index=False)
pd.concat([d for _, b in res for d in b], ignore_index=True).to_csv(f'{R}/ecdf_by_group{SUF}.tsv', sep='\t', index=False)
log('KS rows', len(KS))
print(KS[KS.group == 'WB'][['feature', 'n_esnp', 'n_non', 'median_esnp', 'median_non', 'D', 'direction', 'p']].sort_values('D', ascending=False).to_string())
cx = KS[KS.level == 'ctx']
print('\nctx-level: tested', cx.groupby('feature').size().max(), 'contexts | q<0.05 per feature:')
print(cx.assign(sig=cx.q < 0.05).groupby('feature').agg(n=('sig', 'size'), sig=('sig', 'sum'), median_D=('D', 'median')).sort_values('median_D', ascending=False).to_string())

if ESET != 'eg':
    log('S19_DONE'); sys.exit(0)
pq = pd.read_csv(f'{O}/sets/esnp_pairs_q.tsv.gz', sep='\t', usecols=['ctx', 'dist', 'qval'])
pq = pq[pq.qval < 0.05]
pq['ad'] = pq.dist.abs()
tss = pq.groupby('ctx').agg(n_pairs=('ad', 'size'), pct_10kb=('ad', lambda v: 100 * (v <= 10000).mean()),
                            pct_100kb=('ad', lambda v: 100 * (v <= 100000).mean()), median_abs_dist=('ad', 'median')).reset_index()
tss.to_csv(f'{R}/tss_by_ctx.tsv', sep='\t', index=False)
log('tss_by_ctx', len(tss), 'contexts | pct_100kb range', round(tss.pct_100kb.min(), 1), '-', round(tss.pct_100kb.max(), 1))
log('S19_DONE')
