# -*- coding: utf-8 -*-
"""Why is the sign concordance with Bryois below 50%? Check allele labels on both sides against 1000G EUR.
Run after: plink --bfile g1000_eur --freq --out $OUT/g1000  (A1 = minor allele, MAF = its frequency)."""
import glob
import os

import numpy as np
import pandas as pd

IN = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data'
BRY = '/public/home/yucong/scbrain/external/bryois2022'
OUT = '/public/home/yucong/scbrain/sn_replication'
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def comp(s):
    return s.map(COMP)


g = pd.read_csv(OUT + '/g1000.frq', sep=r'\s+', usecols=['SNP', 'A1', 'A2', 'MAF'])
g.columns = ['snp', 'g1', 'g2', 'gf']
g = g.drop_duplicates('snp', keep='last')
bim = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[1, 4, 5], names=['snp', 'a1', 'a2']).drop_duplicates('snp', keep='last')
frq = pd.read_csv(IN + '/genotype/MAF.txt', sep=r'\s+', usecols=[1, 3, 5], names=['snp', 'alt', 'f1'], header=0).drop_duplicates('snp', keep='last')
pos = pd.read_csv(BRY + '/snp_pos.txt', sep='\t', usecols=['SNP', 'effect_allele', 'other_allele'])
pos.columns = ['snp', 'ea', 'oa']
pos = pos.drop_duplicates('snp', keep='last')


def classify(d, x1, x2, y1, y2):
    same = (d[x1] == d[y1]) & (d[x2] == d[y2])
    swap = (d[x1] == d[y2]) & (d[x2] == d[y1])
    c1, c2 = comp(d[x1]), comp(d[x2])
    fs = (c1 == d[y1]) & (c2 == d[y2]) & ~same & ~swap
    fw = (c1 == d[y2]) & (c2 == d[y1]) & ~same & ~swap
    pal = (d[x1] + d[x2]).isin(['AT', 'TA', 'CG', 'GC'])
    return np.select([pal & (same | swap), same, swap, fs, fw], ['palindromic', 'same', 'swapped', 'strand_same', 'strand_swapped'], 'other')


print('=== 1. our genotype file vs 1000G EUR (all variants) ===')
d = bim.merge(frq, on='snp').merge(g, on='snp')
d['cls'] = classify(d, 'a1', 'a2', 'g1', 'g2')
print('variants in both', len(d), '| ALT label equals bim col 5:', int((d.alt == d.a1).sum()))
for k, x in d.groupby('cls'):
    if k in ('same', 'strand_same'):
        exp = x.gf
    elif k in ('swapped', 'strand_swapped'):
        exp = 1 - x.gf
    else:
        exp = None
    msg = '%-15s n=%8d' % (k, len(x))
    if exp is not None and len(x) > 10:
        msg += ' | corr(freq of our col-5 allele, 1000G freq of the same allele) = %.4f | share with |diff|>0.2: %.4f' % (
            np.corrcoef(x.f1, exp)[0, 1], float((np.abs(x.f1 - exp) > 0.2).mean()))
    print(msg)
x = d[d.cls == 'palindromic']
if len(x):
    e_same = np.where(x.a1 == x.g1, x.gf, 1 - x.gf)
    print('palindromic: corr with 1000G assuming same strand = %.4f' % np.corrcoef(x.f1, e_same)[0, 1])

print('\n=== 2. Bryois allele table vs 1000G EUR and vs our file ===')
b = pos.merge(g, on='snp')
b['cls'] = classify(b, 'ea', 'oa', 'g1', 'g2')
print('Bryois SNPs in 1000G', len(b))
print(b.cls.value_counts().to_string())
print('effect_allele is the 1000G minor allele (non-palindromic, same strand): %.4f' % float((b[b.cls.isin(['same', 'swapped'])].cls == 'same').mean()))
ob = bim.merge(frq, on='snp').merge(pos, on='snp')
ob['cls'] = classify(ob, 'a1', 'a2', 'ea', 'oa')
print('our variants in Bryois table', len(ob))
print(ob.cls.value_counts().to_string())
k = ob[ob.cls.isin(['same', 'swapped'])]
print('Bryois effect_allele equals our col-5 allele: %.4f | our col-5 allele is the minor allele: %.4f' % (
    float((k.cls == 'same').mean()), float((k.f1 < 0.5).mean())))
print('Bryois effect_allele is our minor allele: %.4f' % float(np.where(k.cls == 'same', k.f1 < 0.5, k.f1 > 0.5).mean()))

print('\n=== 3. sign concordance of lead pairs, by subset ===')
lead = pd.read_csv(OUT + '/pairs_lead.tsv.gz', sep='\t').merge(pos, on='snp')
lead['cls'] = classify(lead, 'a1', 'a2', 'ea', 'oa')
bry = []
for f in sorted(glob.glob(OUT + '/bry/*.tsv')):
    if os.path.getsize(f):
        t = pd.read_csv(f, sep='\t', header=None, names=['gene', 'snp', 'dist_b', 'p_b', 'beta_b'])
        t['set'] = os.path.basename(f)[:-4].rsplit('.', 1)[0]
        bry.append(t)
bry = pd.concat(bry).drop_duplicates(['set', 'gene', 'snp'])
m = lead.merge(bry, on=['gene', 'snp'])
m = m[m.cls.isin(['same', 'swapped'])].copy()
m['s'] = np.where(m.cls == 'same', m.slope, -m.slope)
m['agree'] = np.sign(m.s) == np.sign(m.beta_b)


def show(tag, x):
    if len(x) >= 20:
        print('%-46s n=%6d  agree=%.3f' % (tag, len(x), float(x.agree.mean())))


for st in ['pb', 'Excitatory.neurons', 'Oligodendrocytes', 'Astrocytes']:
    x = m[m.set == st]
    print('--', st, '| all pairs n=%d' % len(x))
    show('P_bryois<0.05', x[x.p_b < 0.05])
    show('P_bryois<1e-5', x[x.p_b < 1e-5])
    show('P_bryois<1e-5 and our P<1e-8', x[(x.p_b < 1e-5) & (x.p < 1e-8)])
    show('P_bryois<0.05, labels same', x[(x.p_b < 0.05) & (x.cls == 'same')])
    show('P_bryois<0.05, labels swapped', x[(x.p_b < 0.05) & (x.cls == 'swapped')])
    show('P_bryois<0.05, our col-5 allele minor', x[(x.p_b < 0.05) & (x.af < 0.5)])
    show('P_bryois<0.05, our col-5 allele major', x[(x.p_b < 0.05) & (x.af >= 0.5)])
    show('P_bryois>0.5 (null, expect 0.5)', x[x.p_b > 0.5])
# distance sanity: Bryois distance vs ours for the same pair
x = m[m.set == 'pb']
print('\ncorr(our TSS distance, Bryois TSS distance) = %.4f ; sign agreement of distances = %.4f' % (
    np.corrcoef(x.dist, x.dist_b)[0, 1], float((np.sign(x.dist) == np.sign(x.dist_b)).mean())))
print('DIAG_DONE')
