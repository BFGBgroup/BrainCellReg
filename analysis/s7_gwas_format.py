"""新换的 GWAS 汇总数据转成三种输入，口径照贺奕捷的 3-new.sh / 5-mk-coloc.sh / 7-coloc-non-NA.sh：
  <编号>-smr.tsv      SNP A1 A2 freq b se p n（n 填 NA；频率缺失填 0.5，与他的做法一致，SMR 用 --diff-freq 1）
  <编号>-coloc.tsv    同一张表按 P 升序、每个 SNP 留 P 最小的一行
  <编号>.sumstats.txt SNP P 两列（scGWAS / MAGMA 用）
被删的行逐类计数写进 qc/<编号>.json，不静默丢。用法: python s7_gwas_format.py <编号>
"""
import os, sys, gzip, json, time
import numpy as np
import pandas as pd
from scipy.stats import norm

SRC = '/public/home/yucong/scbrain/external/gwas10'
OUT = os.environ.get('GWAS_OUT', '/public/home/yucong/scbrain/gwas10')
ESI = '/public/home/yucong/scbrain_in/SMR_inputs/SMR/CRBL.Astro.esi'
BIM = '/public/home/yucong/scbrain/local_handoff/scbrain1/g1000_eur/g1000_eur.bim'
SPEC = {
    'GCST90432156': dict(name='MG', path='GCST90432156_MG/GCST90432156.h.tsv.gz', kind='ssf', eff='odds_ratio'),
    # 多发性硬化：Catalog 的 harmonised 文件不能用（元数据把坐标标成 0-based，流程给每个位置加了 1，
    # 只有 21.85% 的位点"对上"且对到的是邻位的 rsID，和 eQTL 的 SNP 只重合 0.56%）。
    # 改用作者原始文件（GRCh37、无 rsID），按 染色体:位置 + 等位基因 到 1000G 欧洲参考里取 rsID。
    'GCST90558092': dict(name='MS', path='GCST90558092_MS/GCST90558092.tsv', kind='ssf_raw_b37'),
    # 10-10 用户定：多发性硬化改用德国队列（Andlauer 2016，4,888 例 / 10,395 对照）；上面那份不用，留着备查
    'GCST003566': dict(name='MS_German', path='GCST003566_MS_German/GCST003566.h.tsv.gz', kind='ssf', eff='beta'),
    'GCST90726344': dict(name='MDD', path='GCST90726344_MDD/GCST90726344.h.tsv.gz', kind='ssf', eff='beta'),
    'GCST009325': dict(name='PD', path='GCST009325_PD/GCST009325.h.tsv.gz', kind='ssf', eff='beta'),
    'GCST90018844': dict(name='FebrileSeizures', path='GCST90018844_FebrileSeizures/GCST90018844.h.tsv.gz', kind='ssf', eff='beta'),
    # PGC 的两份不在 GWAS Catalog 的下载区，编号用 Catalog 里对应研究的编号做文件名
    'GCST90128471': dict(name='SCZ_PGC3_EUR', path='PGC3_SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz', kind='pgc_scz'),
    'GCST011179': dict(name='CUD_PGC_EUR_casecontrol', path='PGC_CUD/CUD_EUR_casecontrol_public_11.14.2020.gz', kind='pgc_cud'),
}
NA = ['NA', '#NA', 'NaN', 'nan', '', '.']
gid = sys.argv[1]
sp = SPEC[gid]
src = os.path.join(SRC, sp['path'])
os.makedirs(os.path.join(OUT, 'qc'), exist_ok=True)
t0 = time.time()
qc = dict(id=gid, name=sp['name'], source=sp['path'])

if sp['kind'] == 'ssf':
    eff = sp['eff']
    cols = ['rsid', 'effect_allele', 'other_allele', 'effect_allele_frequency', eff, 'standard_error', 'p_value']
    d = pd.read_csv(src, sep='\t', usecols=cols, na_values=NA, keep_default_na=False, low_memory=False,
                    dtype={'rsid': str, 'effect_allele': str, 'other_allele': str})
    d = d.rename(columns={'rsid': 'SNP', 'effect_allele': 'A1', 'other_allele': 'A2', 'effect_allele_frequency': 'freq',
                          'standard_error': 'se', 'p_value': 'p'})
    if eff == 'odds_ratio':
        orv = pd.to_numeric(d[eff], errors='coerce')
        d['b'] = np.where(orv > 0, np.log(orv.where(orv > 0)), np.nan)
        qc['effect'] = 'ln(odds_ratio)'
    else:
        d['b'] = pd.to_numeric(d[eff], errors='coerce')
        qc['effect'] = 'beta'
elif sp['kind'] == 'ssf_raw_b37':
    d = pd.read_csv(src, sep='\t', na_values=NA, keep_default_na=False, low_memory=False,
                    dtype={'chromosome': str, 'effect_allele': str, 'other_allele': str})
    qc['rows_raw'] = int(len(d))
    bim = pd.read_csv(BIM, sep='\t', header=None, usecols=[0, 1, 3, 4, 5], names=['chromosome', 'SNP', 'pos', 'ba1', 'ba2'],
                      dtype={'chromosome': str, 'SNP': str, 'ba1': str, 'ba2': str})
    d['ea'] = d.effect_allele.str.upper(); d['oa'] = d.other_allele.str.upper()
    best = None
    for off in (0, 1):      # 两种坐标假设都试，按等位基因对得上的行数定
        x = d.assign(pos=d.base_pair_location + off).merge(bim, on=['chromosome', 'pos'], how='inner')
        okm = ((x.ea == x.ba1) & (x.oa == x.ba2)) | ((x.ea == x.ba2) & (x.oa == x.ba1))
        qc[f'raw_pos_match_offset{off}'] = int(len(x)); qc[f'raw_allele_match_offset{off}'] = int(okm.sum())
        if best is None or okm.sum() > best[0]:
            best = (int(okm.sum()), off, x[okm])
    qc['raw_offset_used'] = best[1]
    qc['drop_not_in_1000G_EUR_or_allele_mismatch'] = int(qc['rows_raw'] - best[2][['chromosome', 'base_pair_location', 'ea', 'oa']].drop_duplicates().shape[0])
    d = best[2].rename(columns={'ea': 'A1', 'oa': 'A2', 'effect_allele_frequency': 'freq', 'beta': 'b', 'standard_error': 'se', 'p_value': 'p'})
    del bim, best
    qc['effect'] = 'beta'
elif sp['kind'] == 'pgc_scz':
    skip = 0
    with gzip.open(src, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip += 1
            else:
                break
    d = pd.read_csv(src, sep='\t', skiprows=skip, na_values=NA, keep_default_na=False, low_memory=False,
                    usecols=['ID', 'A1', 'A2', 'FCAS', 'FCON', 'BETA', 'SE', 'PVAL', 'NCAS', 'NCON'], dtype={'ID': str, 'A1': str, 'A2': str})
    for c in ['FCAS', 'FCON', 'NCAS', 'NCON']:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d['freq'] = (d.FCAS * d.NCAS + d.FCON * d.NCON) / (d.NCAS + d.NCON)
    qc['max_ncas'] = float(d.NCAS.max()); qc['max_ncon'] = float(d.NCON.max())
    d = d.rename(columns={'ID': 'SNP', 'BETA': 'b', 'SE': 'se', 'PVAL': 'p'})
    qc['effect'] = 'BETA = ln(OR) of A1'
elif sp['kind'] == 'pgc_cud':
    d = pd.read_csv(src, sep=r'\s+', na_values=NA, keep_default_na=False, low_memory=False, dtype={'SNP': str, 'A1': str, 'A2': str})
    qc['max_ncas'] = float(pd.to_numeric(d.N_CAS, errors='coerce').max()); qc['max_ncon'] = float(pd.to_numeric(d.N_CON, errors='coerce').max())
    d = d.rename(columns={'Beta': 'b', 'SE': 'se', 'P': 'p'})
    d['freq'] = np.nan          # 公开版不给频率
    qc['effect'] = 'Beta = ln(OR) of A1'
else:
    raise SystemExit('unknown kind')

for c in ['freq', 'b', 'se', 'p']:
    d[c] = pd.to_numeric(d[c], errors='coerce')
d = d[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p']]
qc['rows_in'] = int(len(d))

m1 = d.SNP.isna() | ~d.SNP.astype(str).str.startswith('rs')
qc['drop_no_rsid'] = int(m1.sum()); d = d[~m1]
m2 = d.A1.isna() | d.A2.isna()
qc['drop_allele_missing'] = int(m2.sum()); d = d[~m2]
m3 = ~np.isfinite(d.b) | ~np.isfinite(d.se) | ~np.isfinite(d.p)
qc['drop_b_se_p_missing'] = int(m3.sum()); d = d[~m3]
m4 = d.se <= 0
qc['drop_se_nonpositive'] = int(m4.sum()); d = d[~m4]
m5 = (d.p < 0) | (d.p > 1)
qc['drop_p_out_of_range'] = int(m5.sum()); d = d[~m5]
z0 = d.p == 0
qc['p_zero_floored_to_1e-300'] = int(z0.sum())
d.loc[z0, 'p'] = 1e-300
fm = d.freq.isna()
qc['freq_missing_set_0.5'] = int(fm.sum())
d.loc[fm, 'freq'] = 0.5
d['A1'] = d.A1.str.upper(); d['A2'] = d.A2.str.upper()
d['n'] = 'NA'
qc['rows_smr'] = int(len(d))

# P 与 b/se 是否自洽（判断 se 是不是 b 所在尺度上的标准误）
pc = 2 * norm.sf(np.abs(d.b / d.se))
ok = (d.p > 1e-290) & (pc > 0)
dl = np.abs(np.log10(d.p[ok]) - np.log10(pc[ok]))
qc['p_vs_bse_median_abs_dlog10'] = float(np.median(dl)); qc['p_vs_bse_frac_dlog10_gt0.3'] = float((dl > 0.3).mean())

d.to_csv(os.path.join(OUT, f'{gid}-smr.tsv'), sep='\t', index=False, na_rep='NA', float_format='%.8g')
c = d.sort_values('p', kind='mergesort').drop_duplicates('SNP', keep='first')
qc['rows_coloc'] = int(len(c)); qc['dup_snp_rows_removed'] = int(len(d) - len(c))
qc['n_p_lt_5e-8'] = int((c.p < 5e-8).sum()); qc['min_p'] = float(c.p.min())
qc['median_se'] = float(c.se.median()); qc['median_abs_b'] = float(c.b.abs().median())
c.to_csv(os.path.join(OUT, f'{gid}-coloc.tsv'), sep='\t', index=False, na_rep='NA', float_format='%.8g')
c[['SNP', 'p']].rename(columns={'p': 'P'}).to_csv(os.path.join(OUT, f'{gid}.sumstats.txt'), sep='\t', index=False, float_format='%.8g')

# 与 eQTL 的 SNP 有多少对得上（拿一个组合的 esi 当代表，各组合 SNP 集合相同）
esi = set(pd.read_csv(ESI, sep='\t', header=None, usecols=[1], dtype=str)[1])
qc['eqtl_snps'] = len(esi); qc['eqtl_snps_in_gwas'] = int(c.SNP.isin(esi).sum())
qc['eqtl_snps_in_gwas_frac'] = round(qc['eqtl_snps_in_gwas'] / len(esi), 4)
qc['seconds'] = round(time.time() - t0)
json.dump(qc, open(os.path.join(OUT, 'qc', f'{gid}.json'), 'w'), indent=1)
print(json.dumps(qc, indent=1))
print('FORMAT_DONE', gid)
