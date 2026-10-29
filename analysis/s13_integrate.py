"""把 SMR、共定位、scGWAS 三路结果按「病 × 基因」并到一张表（对应交接文档 1.1 的第二问：
哪些基因在疾病 GWAS 与细胞类型 eQTL 之间同时有 SMR 关联、共定位证据或 scGWAS 模块证据）。
只做汇总计数，不加新的筛选规则：SMR 用 q_SMR<0.05（他的口径）；HEIDI、PP.H4 的各档并列给出。
"""
import pandas as pd

B = '/public/home/yucong/scbrain'
GTF = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data/expression/Homo_sapiens.GRCh38.82.chr.gtf.txt'
NAME = {'GCST90432156': 'MG', 'GCST003566': 'MS', 'GCST90726344': 'MDD', 'GCST009325': 'PD', 'GCST90018844': 'FebrileSeizures',
        'GCST90128471': 'SCZ', 'GCST011179': 'CUD', 'GCST013196': 'AD', 'GCST90014122': 'LacunarStroke', 'GCST90477562': 'Narcolepsy'}
ORDER = ['MG', 'MS', 'AD', 'SCZ', 'MDD', 'LacunarStroke', 'PD', 'FebrileSeizures', 'Narcolepsy', 'CUD']
import os
# 默认读原来的目录；用改正等位基因后的 SMR 结果时由环境变量指到 smr10_fix / coloc10_fix / integrate10_fix
SMR_DIR = os.environ.get('SMR_DIR', f'{B}/smr10')
COLOC_DIR = os.environ.get('COLOC_DIR', f'{B}/coloc10')
OUT = os.environ.get('INT_OUT', f'{B}/integrate10')
os.makedirs(OUT, exist_ok=True)
print('SMR', SMR_DIR, '| coloc', COLOC_DIR, '| out', OUT)

g = pd.read_csv(GTF, sep='\t', usecols=[0, 1, 2, 4, 5, 9], header=0, names=['chr', 'start', 'end', 'gene', 'symbol', 'biotype'], dtype=str)
g = g.drop_duplicates('gene').set_index('gene')
smr = pd.read_csv(f'{SMR_DIR}/smr_q05.tsv', sep='\t')
col = pd.read_csv(f'{COLOC_DIR}/coloc_all.tsv', sep='\t').rename(columns={'ctx': 'Region-Cell', 'gwas': 'GWAS', 'gene': 'Gene'})
m = smr.merge(col[['Region-Cell', 'GWAS', 'Gene', 'nsnps', 'PP.H3.abf', 'PP.H4.abf']], on=['Region-Cell', 'GWAS', 'Gene'], how='left')
m['disease'] = m.GWAS.map(NAME)
m['symbol'] = m.Gene.map(g.symbol); m['biotype'] = m.Gene.map(g.biotype)
m['gene_chr'] = m.Gene.map(g.chr); m['gene_start'] = pd.to_numeric(m.Gene.map(g.start), errors='coerce')
m['MHC'] = (m.gene_chr == '6') & (m.gene_start >= 25e6) & (m.gene_start <= 35e6)        # hg38
m['region'] = m['Region-Cell'].str.split('-').str[0]; m['cell'] = m['Region-Cell'].str.split('-').str[1]
m['heidi_ok_001'] = m.p_HEIDI >= 0.01; m['heidi_ok_005'] = m.p_HEIDI >= 0.05
m['h4_08'] = m['PP.H4.abf'] >= 0.8; m['h4_05'] = m['PP.H4.abf'] >= 0.5
m.to_csv(f'{OUT}/smr_coloc_triplets.tsv', sep='\t', index=False)

# scGWAS：每个病，哪些基因落在显著模块里（任一组合）
sg = pd.read_csv(f'{B}/scgwas10/scgwas_sig_modules_new.tsv', sep='\t')
sgg = {}
for (gw), d in sg.groupby('GWAS'):
    s = {}
    for _, r in d.iterrows():
        for x in str(r.module_genes).split(':'):
            s.setdefault(x, set()).add(r.ctx)
    sgg[gw] = s

rows = []
for (gw, ge), d in m.groupby(['GWAS', 'Gene']):
    sc = sgg.get(gw, {}).get(ge, set())
    rows.append(dict(disease=NAME[gw], GWAS=gw, gene=ge, symbol=d.symbol.iloc[0], biotype=d.biotype.iloc[0], chr=d.gene_chr.iloc[0],
                     start=d.gene_start.iloc[0], MHC=bool(d.MHC.iloc[0]), n_ctx_smr=len(d), n_regions=d.region.nunique(), n_cells=d.cell.nunique(),
                     n_ctx_heidi_ge_0_01=int(d.heidi_ok_001.sum()), n_ctx_h4_ge_0_8=int(d.h4_08.sum()), n_ctx_h4_ge_0_5=int(d.h4_05.sum()),
                     n_ctx_smr_heidi_h4=int((d.heidi_ok_001 & d.h4_08).sum()), max_h4=d['PP.H4.abf'].max(), min_q_smr=d.q_SMR.min(),
                     cells=';'.join(sorted(d.cell.unique())), regions=';'.join(sorted(d.region.unique())),
                     scgwas_module=len(sc) > 0, scgwas_ctx=';'.join(sorted(sc))))
G = pd.DataFrame(rows)
G['disease'] = pd.Categorical(G.disease, ORDER)
G = G.sort_values(['disease', 'n_ctx_smr_heidi_h4', 'n_ctx_h4_ge_0_8', 'n_ctx_smr'], ascending=[True, False, False, False])
G.to_csv(f'{OUT}/gene_evidence.tsv', sep='\t', index=False)

s = []
for dz in ORDER:
    t = m[m.disease == dz]; x = G[G.disease == dz]
    gw = [k for k, v in NAME.items() if v == dz][0]
    s.append(dict(disease=dz, smr_q05_triplets=len(t), smr_genes=len(x), smr_genes_MHC=int(x.MHC.sum()), smr_genes_nonMHC=int((~x.MHC).sum()),
                  triplets_heidi_ge_0_01=int(t.heidi_ok_001.sum()), triplets_h4_ge_0_8=int(t.h4_08.sum()),
                  genes_h4_ge_0_8=int((x.n_ctx_h4_ge_0_8 > 0).sum()), genes_h4_ge_0_8_nonMHC=int(((x.n_ctx_h4_ge_0_8 > 0) & ~x.MHC).sum()),
                  genes_smr_heidi_h4=int((x.n_ctx_smr_heidi_h4 > 0).sum()), genes_smr_heidi_h4_nonMHC=int(((x.n_ctx_smr_heidi_h4 > 0) & ~x.MHC).sum()),
                  scgwas_ctx_with_modules=sg[sg.GWAS == gw].ctx.nunique(), scgwas_module_genes=len(sgg.get(gw, {})),
                  genes_smr_and_scgwas=int(x.scgwas_module.sum()), genes_all_three=int(((x.n_ctx_h4_ge_0_8 > 0) & x.scgwas_module).sum())))
S = pd.DataFrame(s)
S.to_csv(f'{OUT}/summary_by_disease.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(S.to_string(index=False))
print('\n总计：SMR q<0.05 的（组合、病、基因）', len(m), '；病 × 基因', len(G), '；PP.H4>=0.8 的三元组', int(m.h4_08.sum()),
      '；SMR + HEIDI>=0.01 + PP.H4>=0.8 的三元组', int((m.heidi_ok_001 & m.h4_08).sum()),
      '；对应病 × 基因', int((G.n_ctx_smr_heidi_h4 > 0).sum()), '其中不在 MHC 的', int(((G.n_ctx_smr_heidi_h4 > 0) & ~G.MHC).sum()))
print('\n每种病证据最全的基因（SMR + HEIDI>=0.01 + PP.H4>=0.8 的组合数排前 8）：')
for dz in ORDER:
    x = G[(G.disease == dz) & (G.n_ctx_smr_heidi_h4 > 0)].head(8)
    print(dz, '|', '; '.join(f"{r.symbol}{'*' if r.MHC else ''}({r.n_ctx_smr_heidi_h4}{'+sc' if r.scgwas_module else ''})" for r in x.itertuples()))
print('（* = 基因在 6 号染色体 25–35 Mb；+sc = 同一种病的 scGWAS 显著模块里也有它）')
print('INTEGRATE_DONE')
