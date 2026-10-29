"""Supplementary Data 3 (scGWAS significant modules) and 4 (external SMR replication), column renames and gene symbols
only, no numbers changed. Sources: scgwas10/scgwas_sig_modules_new.tsv (+ scgwas_summary_new.tsv for disease names),
indep/smr_ext_summary.tsv and indep/smr_ext_pairs.tsv. Writes supp/data/.
/public/home/yucong/.conda/envs/scvi/bin/python s82_supp_data.py"""
import gzip, os, re
import duckdb
import pandas as pd

S = '/public/home/yucong/scbrain'
OUT = f'{S}/supp/data'
os.makedirs(OUT, exist_ok=True)

# gene symbols: tested genes first, then the Ensembl 82 gene records for module genes outside the eQTL gene set
sym = dict(duckdb.connect().execute(f"SELECT DISTINCT gene, symbol FROM '{S}/db/tables/egene_stats.parquet'").fetchall())
gtf = f'{S}/local_handoff/scbrain3/Homo_sapiens.GRCh38.82.chr(only_gene).gtf'
with open(gtf, encoding='utf-8', errors='replace') as f:
    for line in f:
        m = re.search(r'gene_id "(ENSG\d+)".*?gene_name "([^"]+)"', line)
        if m and m.group(1) not in sym:
            sym[m.group(1)] = m.group(2)

# ---- Supplementary Data 3
md = pd.read_csv(f'{S}/scgwas10/scgwas_sig_modules_new.tsv', sep='\t')
sm = pd.read_csv(f'{S}/scgwas10/scgwas_summary_new.tsv', sep='\t')
dis = sm[['GWAS', 'disease']].drop_duplicates().set_index('GWAS')['disease']
assert md.GWAS.isin(dis.index).all()
genes = md.module_genes.str.split(':')
miss = sorted({g for gl in genes for g in gl if g not in sym} | {g for g in md.seed if g not in sym})
d3 = pd.DataFrame({
    'context': md.ctx, 'disease': md.GWAS.map(dis), 'gwas_accession': md.GWAS,
    'seed_gene': md.seed, 'seed_symbol': md.seed.map(sym),
    'n_module_genes': genes.str.len(),
    'module_genes_ensembl': md.module_genes.str.replace(':', ';'),
    'module_genes_symbol': genes.apply(lambda gl: ';'.join(sym.get(g, g) for g in gl)),
    'reference_cell_type': md.cell_type,
    'module_score': md.module_score, 'gwas_score': md.m_gwas, 'expression_score': md.m_scrnaseq,
    'module_score_z': md.module_score_z, 'gwas_score_z': md.z_gwas, 'expression_score_z': md.z_scrnaseq,
    'p_module': md.p, 'p_gwas': md.p_gwas, 'p_expression': md.p_scrn})
d3 = d3.sort_values(['disease', 'context', 'p_module']).reset_index(drop=True)
d3.to_csv(f'{OUT}/Supplementary_Data_3_scGWAS_modules.tsv', sep='\t', index=False)
print('Data 3 rows', len(d3), '| combinations', d3[['context', 'disease']].drop_duplicates().shape[0],
      '| diseases', d3.disease.nunique(), sorted(d3.disease.unique()), '| symbols missing', len(miss))

# ---- Supplementary Data 4
su = pd.read_csv(f'{S}/indep/smr_ext_summary.tsv', sep='\t')
pr = pd.read_csv(f'{S}/indep/smr_ext_pairs.tsv', sep='\t')
pr['dataset'] = pr.dataset.str.replace(r'^BrainMeta_cis_eQTL_chr\d+$', 'BrainMeta v2', regex=True)
d4a = su.rename(columns={'gene': 'gene_ensembl', 'all3': 'passed_SMR_HEIDI_coloc', 'datasets_tested': 'n_datasets_tested',
                         'datasets_rep': 'n_datasets_replicated', 'best_p': 'best_external_p_SMR', 'brainmeta_p': 'brainmeta_p_SMR',
                         'any_heidi_ok': 'any_replicating_dataset_p_HEIDI_ge_0.01'})
d4b = pr[['disease', 'gene', 'symbol', 'dataset', 'topSNP', 'A1', 'A2', 'Freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL', 'se_eQTL',
          'p_eQTL', 'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI', 'nsnp_HEIDI', 'b_main', 'p_main', 'same_sign', 'rep']].rename(
    columns={'gene': 'gene_ensembl', 'b_main': 'b_SMR_this_study', 'p_main': 'p_SMR_this_study', 'same_sign': 'same_direction',
             'rep': 'replicated'})
d4a.to_csv(f'{OUT}/Supplementary_Data_4_external_SMR_pairs.tsv', sep='\t', index=False)
d4b.to_csv(f'{OUT}/Supplementary_Data_4_external_SMR_datasets.tsv', sep='\t', index=False)
t = d4a[d4a.n_datasets_tested > 0]
for k, g in t.groupby('passed_SMR_HEIDI_coloc'):
    print('Data 4 all3 =', k, '| pairs', len(d4a[d4a.passed_SMR_HEIDI_coloc == k]), '| testable', len(g), '| replicated', int(g.replicated.sum()))
print('Data 4 dataset rows', len(d4b), '| datasets', d4b.dataset.nunique())
print('DATA_DONE')
