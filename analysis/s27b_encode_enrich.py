"""Fig 4 (figure plan v2) inputs that are not yet in esnp_eg/results:
  1. eSNP enrichment in prefrontal-cortex histone-mark peaks (7 marks) and DNase peaks from ENCODE (reference manuscript Fig 4c/e).
     Peaks of all experiments on donors without a neurological diagnosis are pooled per feature (sort + merge);
     a variant is "in" a feature when its position falls in a merged peak. Fisher test eSNP vs non-eSNP, as step10 (s4_stats.py).
  2. MAF of eSNPs vs non-eSNPs by cell type: ECDF on a fixed grid + two-sided KS (reference manuscript Fig 3d).
Groups and sets exactly as s4_stats.py (SCB_OUT = esnp_eg: eSNPs of eGenes; non-eSNPs = nominal P > 1e-4 everywhere in the group).
Needs bedtools (module apps/bedtools/2.30.0). Output: esnp_eg/results/encode_pfc_enrich.tsv, maf_ecdf_by_group.tsv, maf_ks_by_group.tsv"""
import os, subprocess, glob
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
GROUPS = {'WB': list(range(len(ctx)))}
for c in CELLS:
    GROUPS['cell:' + c] = [j for j, x in enumerate(ctx) if x.split('.')[1] == c]
for r in REGIONS:
    GROUPS['region:' + r] = [j for j, x in enumerate(ctx) if x.split('.')[0] == r]
MK = {}
for g, J in GROUPS.items():
    Eg = e3[:, J].any(1)
    MK[g] = (Eg, n4[:, J].any(1) & ~pre[:, J].any(1) & ~Eg)
print('WB eSNPs', int(MK['WB'][0].sum()), 'non-eSNPs', int(MK['WB'][1].sum()), flush=True)

# ---- 1. ENCODE prefrontal cortex peaks
man = pd.read_csv(f'{E}/manifest.tsv', sep='\t').fillna('')
bad = man.biosample.str.contains('Alzheimer|cognitive impairment|dementia', case=False, regex=True)
use = man[(man.file != '') & ~bad].copy()
use['path'] = [f"{E}/{r.feature}__{r.experiment}__{r.file}.bed.gz" for r in use.itertuples()]
use = use[use.path.map(os.path.exists)]
print('experiments in manifest', len(man), '| excluded for diagnosis', int(bad.sum()), '| used', len(use), flush=True)
print(use.groupby('feature').size().to_string(), flush=True)
use.to_csv(f'{R}/encode_pfc_used.tsv', sep='\t', index=False)
rs_index = pd.Index(var.rsid)
FLAGS = {}
for feat, sub in use.groupby('feature'):
    merged = f'{E}/{feat}.merged.bed'
    cmd = (f"zcat {' '.join(sub.path)} | cut -f1-3 | grep -E '^chr([0-9]+|X|Y)\\b' | sort -k1,1 -k2,2n | bedtools merge > {merged} && "
           f"bedtools intersect -a {O}/sets/universe_hg38.bed -b {merged} -u | cut -f4")
    hits = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout.split()
    f = np.zeros(len(var), bool); i = rs_index.get_indexer(hits); f[i[i >= 0]] = True
    FLAGS[feat] = f
    print(feat, 'experiments', len(sub), 'variants in peaks', int(f.sum()), flush=True)
out = []
for g, (Eg, N) in MK.items():
    nE, nN = int(Eg.sum()), int(N.sum())
    for feat, f in FLAGS.items():
        a, c = int((Eg & f).sum()), int((N & f).sum())
        orr, lo, hi, p = fisher(a, nE - a, c, nN - c)
        out.append(dict(group=g, feature=feat, n_esnp=nE, n_non=nN, esnp_in=a, non_in=c, OR=orr, OR_lo=lo, OR_hi=hi, p=p))
res = pd.DataFrame(out)
res.to_csv(f'{R}/encode_pfc_enrich.tsv', sep='\t', index=False)
print(res[res.group == 'WB'].to_string(), flush=True)

# ---- 2. MAF by group
maf = var['MAF'].values.astype(float)
grid = np.linspace(0, 0.5, 101)
er, kr = [], []
for g, (Eg, N) in MK.items():
    x = maf[Eg]; y = maf[N]; x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
    if len(x) < 5:
        continue
    for s, v in (('eSNP', x), ('non-eSNP', y)):
        F = np.searchsorted(np.sort(v), grid, side='right') / len(v)
        er.append(pd.DataFrame(dict(group=g, set=s, maf=grid, F=F)))
    ks = stats.ks_2samp(x, y)
    kr.append(dict(group=g, n_esnp=len(x), n_non=len(y), median_esnp=np.median(x), median_non=np.median(y), D=ks.statistic, p=ks.pvalue))
pd.concat(er).to_csv(f'{R}/maf_ecdf_by_group.tsv', sep='\t', index=False)
pd.DataFrame(kr).to_csv(f'{R}/maf_ks_by_group.tsv', sep='\t', index=False)
print(pd.DataFrame(kr).to_string(), flush=True)
print('ENCODE_ENRICH_DONE')
