"""SMR 结果汇总，口径照贺奕捷 3-qval.R / 4-final.sh / 5-all-sig-new.sh：
每个结果文件（一个组合 × 一份 GWAS）内对 p_SMR 做 BH → q_SMR；p_SMR<0.05 的汇成 all_sig_with_origin.txt；
q_SMR<0.05 的（组合、GWAS、基因）写 smr-sig-gene.txt，供共定位用。
新换的 7 份 GWAS 读 smr10/ 下新跑的结果；没换的 3 份读他原来的结果，并和他的 smr-sig-gene.txt 对一遍。
"""
import os, sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

NEW = '/public/home/yucong/scbrain/smr10'
HIS = '/public/home/yucong/scbrain_in/paper2/eqtl-gwas/SMR/10disease-SMR'
NEWIDS = ['GCST90432156', 'GCST003566', 'GCST90726344', 'GCST009325', 'GCST90018844', 'GCST90128471', 'GCST011179']
KEEP = ['GCST013196', 'GCST90014122', 'GCST90477562']
# 暂不纳入的编号（逗号分隔）：多发性硬化 GCST90558092 的文件与 Catalog 记载不符（见 TASKS），数据集定下来之前先不进汇总
EXCL = [x for x in os.environ.get('SMR_EXCLUDE', '').split(',') if x]
NEWIDS = [g for g in NEWIDS if g not in EXCL]
print('excluded', EXCL)
NAME = {'GCST90432156': 'MG', 'GCST003566': 'MS', 'GCST90726344': 'MDD', 'GCST009325': 'PD', 'GCST90018844': 'FebrileSeizures',
        'GCST90128471': 'SCZ', 'GCST011179': 'CUD', 'GCST013196': 'AD', 'GCST90014122': 'LacunarStroke', 'GCST90477562': 'Narcolepsy'}
REG = ["MEDU", "PUTM", "TCTX", "WHMT", "CRBL", "OCTX", "SNIG", "HIPP", "THAL", "FCTX"]
CELL = ["Endo", "Purk", "Perc", "NPC", "OPC", "InN", "ExN", "Gran", "Olig", "Astro", "Micro"]

rows, missing, empty = [], [], []
for r in REG:
    for c in CELL:
        d = f'{r}-{c}'
        for g in NEWIDS + KEEP:
            f = f'{NEW}/{d}/{g}.smr' if g in NEWIDS else f'{HIS}/{d}/{g}.smr'
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
cols = ['Region-Cell', 'GWAS', 'Gene', 'topSNP', 'p_GWAS', 'p_eQTL', 'p_SMR', 'p_HEIDI', 'q_SMR']
a[a.p_SMR < 0.05][cols].to_csv(f'{NEW}/all_sig_with_origin.txt', sep='\t', index=False)
sig = a[a.q_SMR < 0.05]
sig[['Region-Cell', 'GWAS', 'Gene']].to_csv(f'{NEW}/smr-sig-gene.txt', sep='\t', index=False, header=False)
sig[cols + ['b_SMR', 'se_SMR', 'nsnp_HEIDI', 'topSNP_chr', 'topSNP_bp']].to_csv(f'{NEW}/smr_q05.tsv', sep='\t', index=False)

print('files read', len(rows), 'missing', len(missing), missing[:6], 'empty', len(empty))
print('rows', len(a), '| p_SMR<0.05', int((a.p_SMR < 0.05).sum()), '| q_SMR<0.05', len(sig))
s = []
for g in NEWIDS + KEEP:
    x = a[a.GWAS == g]; y = x[x.q_SMR < 0.05]
    s.append(dict(GWAS=g, disease=NAME[g], contexts=x['Region-Cell'].nunique(), tests=len(x), p05=int((x.p_SMR < 0.05).sum()), q05=len(y),
                  q05_heidi_ge_0_05=int((y.p_HEIDI >= 0.05).sum()), q05_heidi_ge_0_01=int((y.p_HEIDI >= 0.01).sum()), q05_heidi_na=int(y.p_HEIDI.isna().sum()),
                  q05_genes=y.Gene.nunique(), q05_contexts=y['Region-Cell'].nunique()))
s = pd.DataFrame(s)
s.to_csv(f'{NEW}/smr_summary_by_gwas.tsv', sep='\t', index=False)
print(s.to_string(index=False))

# 没换的 3 份：和他的 smr-sig-gene.txt 对一遍
his = pd.read_csv(f'{HIS}/smr-sig-gene.txt', sep='\t', header=None, names=['Region-Cell', 'GWAS', 'Gene'])
hk = set(map(tuple, his[his.GWAS.isin(KEEP)].values)); mk = set(map(tuple, sig[sig.GWAS.isin(KEEP)][['Region-Cell', 'GWAS', 'Gene']].values))
print('对回原结果（没换的 3 份，q_SMR<0.05）：他', len(hk), '我', len(mk), '相同', len(hk & mk))
print('SMR_POST_DONE')
