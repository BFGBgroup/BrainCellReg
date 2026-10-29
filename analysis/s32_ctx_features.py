"""Fig 4g atlas of the 110 cell type-brain region combinations (PI #33/#36): the per-combination tests that were not yet
computed. Same eSNP / non-eSNP sets as s4_stats.py and s27b_encode_enrich.py (esnp_eg: eSNPs of eGenes; non-eSNPs =
nominal P > 1e-4 everywhere in the group). Groups: WB, cell:*, region:*, ctx:* (110).
  1. ENCODE prefrontal-cortex histone marks (7) and DNase peaks: Fisher, eSNP vs non-eSNP (merged peaks from s27b)
  2. ReMap SIN3A / OLIG2 peaks in brain-prefrontal-cortex: Fisher (as step10_remap_brain_biotypes in s4_stats.py)
  3. MAF: two-sided KS, eSNP vs non-eSNP (as maf_ks_by_group in s27b)
  4. TSS: eSNP-eGene pairs within 10 kb of the TSS vs all tested pairs of the group (as the 10 kb window of step6_tss_by_strength)
Output: esnp_eg/results/ctx_features_extra.tsv (group, feature, test, n1, n2, a, c, OR, OR_lo, OR_hi, D, direction, p)"""
import os, subprocess
import numpy as np
import pandas as pd
from scipy import stats

S = '/public/home/yucong/scbrain'
O = f'{S}/esnp_eg'
R = f'{O}/results'
E = f'{S}/external/encode_pfc'
REGIONS = ['CRBL', 'FCTX', 'HIPP', 'MEDU', 'OCTX', 'PUTM', 'SNIG', 'TCTX', 'THAL', 'WHMT']
CELLS = ['Astro', 'Endo', 'ExN', 'Gran', 'InN', 'Micro', 'NPC', 'Olig', 'OPC', 'Perc', 'Purk']


def fisher(a, b, c, d):
    a, b, c, d = int(a), int(b), int(c), int(d)
    if (a + b) == 0 or (c + d) == 0 or (a + c) == 0:
        return np.nan, np.nan, np.nan, np.nan
    _, p = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (a + .5, b + .5, c + .5, d + .5) if min(a, b, c, d) == 0 else (a, b, c, d)
    lor = np.log(aa * dd / (bb * cc)); se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se), p


var = pd.read_csv(f'{O}/sets/variants.tsv.gz', sep='\t', dtype={'chr': str})
m = np.load(f'{O}/sets/membership.npz', allow_pickle=True)
ctx = [str(x) for x in m['ctx']]
e3, n4, pre = m['e3'], m['n4'], m['pre']
G = {'WB': list(range(len(ctx)))}
for c in CELLS:
    G['cell:' + c] = [j for j, x in enumerate(ctx) if x.split('.')[1] == c]
for r in REGIONS:
    G['region:' + r] = [j for j, x in enumerate(ctx) if x.split('.')[0] == r]
for j, x in enumerate(ctx):
    G['ctx:' + x] = [j]
MK = {}
for g, J in G.items():
    Eg = e3[:, J].any(1)
    MK[g] = (Eg, n4[:, J].any(1) & ~pre[:, J].any(1) & ~Eg)
print('groups', len(G), '| WB eSNPs', int(MK['WB'][0].sum()), flush=True)
rs_index = pd.Index(var.rsid)
out = []


def fisher_flags(flags, test):
    for g, (Eg, N) in MK.items():
        nE, nN = int(Eg.sum()), int(N.sum())
        for feat, f in flags.items():
            a, c = int((Eg & f).sum()), int((N & f).sum())
            orr, lo, hi, p = fisher(a, nE - a, c, nN - c)
            out.append(dict(group=g, feature=feat, test=test, n1=nE, n2=nN, a=a, c=c, OR=orr, OR_lo=lo, OR_hi=hi, p=p))


# 1. ENCODE (merged peak files written by s27b_encode_enrich.py)
FL = {}
for feat in ['H3K27me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K36me3', 'DNase']:
    mb = f'{E}/{feat}.merged.bed'
    hits = subprocess.run(f"bedtools intersect -a {O}/sets/universe_hg38.bed -b {mb} -u | cut -f4", shell=True,
                          capture_output=True, text=True, check=True).stdout.split()
    f = np.zeros(len(var), bool); i = rs_index.get_indexer(hits); f[i[i >= 0]] = True
    FL[feat] = f
    print(feat, int(f.sum()), flush=True)
fisher_flags(FL, 'fisher')

# 2. ReMap brain-prefrontal-cortex SIN3A / OLIG2
TF = {}
for ch in pd.read_csv(f'{O}/annot/remap_brain_hits.tsv', sep='\t', header=None, names=['rsid', 'tf', 'biotype'], chunksize=2_000_000):
    ch = ch[ch.tf.isin(['SIN3A', 'OLIG2'])]
    ch = ch[ch.biotype.str.split(',').map(lambda b: 'brain-prefrontal-cortex' in b)]
    for t, sub in ch.groupby('tf'):
        TF.setdefault(t + ':brain-prefrontal-cortex', []).extend(sub.rsid.tolist())
FT = {}
for k, rs in TF.items():
    f = np.zeros(len(var), bool); i = rs_index.get_indexer(pd.unique(np.array(rs))); f[i[i >= 0]] = True
    FT[k] = f
    print(k, int(f.sum()), flush=True)
fisher_flags(FT, 'fisher')

# 3. MAF KS
maf = var['MAF'].values.astype(float)
for g, (Eg, N) in MK.items():
    x = maf[Eg]; y = maf[N]; x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
    if len(x) < 5:
        out.append(dict(group=g, feature='MAF', test='ks', n1=len(x), n2=len(y)))
        continue
    ks = stats.ks_2samp(x, y)
    out.append(dict(group=g, feature='MAF', test='ks', n1=len(x), n2=len(y), D=ks.statistic, p=ks.pvalue,
                    direction='higher' if np.median(x) > np.median(y) else 'lower'))

# 4. TSS: eSNP-eGene pairs within 10 kb (bins -1, 0) vs all tested pairs
tested = pd.read_csv(f'{O}/tss/tested_pairs_bins.tsv', sep='\t')
pairs = pd.read_csv(f'{O}/sets/esnp_pairs.tsv.gz', sep='\t', usecols=['ctx', 'dist'])
pairs['bin'] = np.floor(pairs['dist'] / 10000.0).astype(int)
tested['inw'] = tested['bin'].between(-1, 0); pairs['inw'] = pairs['bin'].between(-1, 0)
tt = tested.groupby(['ctx', 'inw'])['n'].sum().unstack(fill_value=0)
pp = pairs.groupby(['ctx', 'inw']).size().unstack(fill_value=0)
for g, J in G.items():
    cs = [ctx[j] for j in J]
    t_in, t_out = int(tt.reindex(cs).get(True, pd.Series(dtype=float)).fillna(0).sum()), int(tt.reindex(cs).get(False, pd.Series(dtype=float)).fillna(0).sum())
    e_in, e_out = int(pp.reindex(cs).get(True, pd.Series(dtype=float)).fillna(0).sum()), int(pp.reindex(cs).get(False, pd.Series(dtype=float)).fillna(0).sum())
    orr, lo, hi, p = fisher(e_in, e_out, t_in - e_in, t_out - e_out)
    out.append(dict(group=g, feature='TSS_10kb', test='fisher', n1=e_in + e_out, n2=t_in + t_out, a=e_in, c=t_in,
                    OR=orr, OR_lo=lo, OR_hi=hi, p=p))

res = pd.DataFrame(out)
res.to_csv(f'{R}/ctx_features_extra.tsv', sep='\t', index=False)
print(res[res.group == 'WB'].to_string(), flush=True)
c = res[res.group.str.startswith('ctx:')]
print('ctx rows', len(c), '| non-NA p', int(c.p.notna().sum()), flush=True)
print('CTX_FEATURES_DONE')
