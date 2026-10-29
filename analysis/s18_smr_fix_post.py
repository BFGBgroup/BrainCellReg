# -*- coding: utf-8 -*-
"""SMR rerun with corrected allele labels (smr10_fix): same summary as s10_smr_post.py (BH within each context x GWAS
file, q_SMR < 0.05 -> smr-sig-gene.txt), then a comparison with the run on the handover .esi (smr10 + his own results
for the three GWAS that were not replaced), and the list of triplets that still need coloc.
"""
import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

B = '/public/home/yucong/scbrain'
NEW = B + '/smr10_fix'
OLD = B + '/smr10/smr_all_results.tsv.gz'
COLOC_OLD = B + '/coloc10/coloc_all.tsv'
COLOC_NEW = B + '/coloc10_fix'
NAME = {'GCST90432156': 'MG', 'GCST003566': 'MS', 'GCST90726344': 'MDD', 'GCST009325': 'PD', 'GCST90018844': 'FebrileSeizures',
        'GCST90128471': 'SCZ', 'GCST011179': 'CUD', 'GCST013196': 'AD', 'GCST90014122': 'LacunarStroke', 'GCST90477562': 'Narcolepsy'}
REG = ["MEDU", "PUTM", "TCTX", "WHMT", "CRBL", "OCTX", "SNIG", "HIPP", "THAL", "FCTX"]
CELL = ["Endo", "Purk", "Perc", "NPC", "OPC", "InN", "ExN", "Gran", "Olig", "Astro", "Micro"]
KEY = ['Region-Cell', 'GWAS', 'Gene']

rows, missing, empty = [], [], []
for r in REG:
    for c in CELL:
        d = f'{r}-{c}'
        for g in NAME:
            f = f'{NEW}/{d}/{g}.smr'
            if not os.path.exists(f):
                missing.append((d, g)); continue
            x = pd.read_csv(f, sep='\t')
            if len(x) == 0:
                empty.append((d, g)); continue
            ok = x.p_SMR.notna().values
            q = np.full(len(x), np.nan)
            if ok.any():
                q[ok] = multipletests(x.p_SMR.values[ok], method='fdr_bh')[1]
            x['q_SMR'] = q
            x.insert(0, 'GWAS', g); x.insert(0, 'Region-Cell', d)
            rows.append(x)
a = pd.concat(rows, ignore_index=True)
a.to_csv(f'{NEW}/smr_all_results.tsv.gz', sep='\t', index=False)
cols = KEY + ['topSNP', 'p_GWAS', 'p_eQTL', 'p_SMR', 'p_HEIDI', 'q_SMR']
a[a.p_SMR < 0.05][cols].to_csv(f'{NEW}/all_sig_with_origin.txt', sep='\t', index=False)
sig = a[a.q_SMR < 0.05]
sig[KEY].to_csv(f'{NEW}/smr-sig-gene.txt', sep='\t', index=False, header=False)
sig[cols + ['b_SMR', 'se_SMR', 'nsnp_HEIDI', 'topSNP_chr', 'topSNP_bp']].to_csv(f'{NEW}/smr_q05.tsv', sep='\t', index=False)
print('files read', len(rows), '| missing', len(missing), missing[:6], '| empty', len(empty))
print('rows', len(a), '| p_SMR<0.05', int((a.p_SMR < 0.05).sum()), '| q_SMR<0.05', len(sig))

o = pd.read_csv(OLD, sep='\t')
osig = o[o.q_SMR < 0.05]
s = []
for g in NAME:
    x, y = a[a.GWAS == g], sig[sig.GWAS == g]
    xo, yo = o[o.GWAS == g], osig[osig.GWAS == g]
    both = y.merge(yo[KEY], on=KEY)
    s.append(dict(GWAS=g, disease=NAME[g], tests_old=len(xo), tests_new=len(x), q05_old=len(yo), q05_new=len(y), q05_both=len(both),
                  genes_old=yo.Gene.nunique(), genes_new=y.Gene.nunique(), genes_both=len(set(yo.Gene) & set(y.Gene)),
                  heidi_ge_0_01_old=int((yo.p_HEIDI >= 0.01).sum()), heidi_ge_0_01_new=int((y.p_HEIDI >= 0.01).sum()),
                  heidi_na_new=int(y.p_HEIDI.isna().sum())))
S = pd.DataFrame(s)
S.to_csv(f'{NEW}/smr_summary_old_vs_new.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(S.to_string(index=False))
print('total q_SMR<0.05: old %d, new %d, both %d | genes (disease x gene) old %d, new %d, both %d' % (
    len(osig), len(sig), int(S.q05_both.sum()), int(S.genes_old.sum()), int(S.genes_new.sum()), int(S.genes_both.sum())))

m = o.merge(a, on=KEY, suffixes=('_o', '_n'))
st = m[m.topSNP_o == m.topSNP_n]
print('tests in both runs %d | same top SNP %.4f | among those: corr of -log10 p_SMR %.5f, b_SMR same sign %.4f' % (
    len(m), len(st) / len(m), np.corrcoef(-np.log10(st.p_SMR_o), -np.log10(st.p_SMR_n))[0, 1],
    float((np.sign(st.b_SMR_o) == np.sign(st.b_SMR_n)).mean())))
ms = m[(m.q_SMR_o < 0.05) | (m.q_SMR_n < 0.05)]
hh = ms.dropna(subset=['p_HEIDI_o', 'p_HEIDI_n'])
print('q_SMR<0.05 in either run and present in both: %d (same top SNP %.4f) | HEIDI in both %d: pass (>=0.01) old %d, new %d, verdict agrees %.4f' % (
    len(ms), float((ms.topSNP_o == ms.topSNP_n).mean()), len(hh), int((hh.p_HEIDI_o >= 0.01).sum()), int((hh.p_HEIDI_n >= 0.01).sum()),
    float(((hh.p_HEIDI_o >= 0.01) == (hh.p_HEIDI_n >= 0.01)).mean())))

# triplets that have no coloc result yet (coloc does not use allele labels, so existing results are reused)
os.makedirs(COLOC_NEW, exist_ok=True)
co = pd.read_csv(COLOC_OLD, sep='\t').rename(columns={'ctx': 'Region-Cell', 'gwas': 'GWAS', 'gene': 'Gene'})
todo = sig[KEY].merge(co[KEY].drop_duplicates(), on=KEY, how='left', indicator=True)
todo = todo[todo._merge == 'left_only'][KEY]
todo.to_csv(f'{COLOC_NEW}/smr-sig-gene.txt', sep='\t', index=False, header=False)
print('coloc: triplets with q_SMR<0.05 %d, already have a result %d, to run %d' % (len(sig), len(sig) - len(todo), len(todo)))
print('SMR_FIX_POST_DONE')
