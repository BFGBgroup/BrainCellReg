"""Database step 1 (PI 10-12: "先生成表格，再加载画图函数……最后生成数据库"): the result tables, one Parquet file each,
under /public/home/yucong/scbrain/db/tables, plus a check of the key numbers against the paper (db/tables/CHECK.txt).
Tables that wait for the open decisions of the database plan are not built here: per-sample genotypes (individual-level
data) and the 176 GB of nominal cis results (queried in place if they go on the server).
  contexts         110 brain region x cell type contexts: samples, mean cell proportion, genes tested, eGenes, independent eSNPs
  samples          1,231 samples: region, individual, sex, age group, 11 estimated cell proportions
  genes            genes tested: Ensembl id, symbol, biotype, position (GRCh38)
  egene_stats      tensorQTL permutation result for every gene in every context (lead variant, effect allele, slope, q)
  esnp_pairs       eSNP-gene pairs of eGenes (q<0.05): position, distance, P, slope, effect allele; independent-eSNP flag
  fisher           Fisher enrichment at GWAS Catalog loci, neurological and non-neurological traits, significance flag
  smr              SMR of the 10 diseases x 110 contexts (corrected alleles), with HEIDI and coloc PP.H4 where tested
  gene_evidence    disease x gene summary (SMR, HEIDI, colocalisation, scGWAS module)
  scgwas_modules   significant scGWAS modules
"""
import os, sys, glob
import numpy as np
import pandas as pd
import duckdb

B = '/public/home/yucong/scbrain'
IN = '/public/home/yucong/scbrain_in/paper2'
OUT = B + '/db/tables'
os.makedirs(OUT, exist_ok=True)
REG = ['FCTX', 'TCTX', 'OCTX', 'HIPP', 'PUTM', 'THAL', 'SNIG', 'MEDU', 'WHMT', 'CRBL']
CELL = ['ExN', 'InN', 'Astro', 'Olig', 'OPC', 'Micro', 'Endo', 'Perc', 'NPC', 'Gran', 'Purk']
con = duckdb.connect()
check = []


def save(df, name):
    con.register('t', df)
    con.execute(f"COPY t TO '{OUT}/{name}.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)")
    con.unregister('t')
    print(f'{name}: {len(df):,} rows, {df.shape[1]} columns', flush=True)


def ok(label, got, want):
    check.append(f"{'OK ' if got == want else 'BAD'} {label}: {got} (paper {want})")


# allele table: counted allele of the genotype dosages (= effect allele of tensorQTL slopes)
at = pd.read_csv(B + '/allele_fix/allele_table_v2.tsv.gz', sep='\t', usecols=['snp', 'eff', 'oth', 'status'])
at = at.drop_duplicates('snp').rename(columns={'snp': 'rsid', 'eff': 'effect_allele', 'oth': 'other_allele', 'status': 'allele_status'})
at.loc[~at.allele_status.isin(['ok', 'ok_ld']), ['effect_allele', 'other_allele']] = None

# genes
gtf = pd.read_csv(IN + '/single-cell-eqtl/data/expression/Homo_sapiens.GRCh38.82.chr.gtf.txt', sep='\t', usecols=[0, 1, 2, 3, 4, 5, 9], header=0,
                  names=['chr', 'start', 'end', 'strand', 'gene', 'symbol', 'biotype'], dtype={'chr': str}).drop_duplicates('gene')

# egene_stats from the permutation pass (run/new2/<ctx>.cis_qtl.txt)
rows = []
for f in sorted(glob.glob(IN + '/single-cell-eqtl/data/run/new2/*.cis_qtl.txt')):
    ctx = os.path.basename(f).replace('.cis_qtl.txt', '')
    d = pd.read_csv(f, sep='\t')
    d.insert(0, 'ctx', ctx)
    rows.append(d)
es = pd.concat(rows, ignore_index=True).rename(columns={'phenotype_id': 'gene', 'variant_id': 'lead_rsid'})
es[['region', 'cell']] = es.ctx.str.split('.', expand=True)
es['is_egene'] = es.qval < 0.05
es = es.merge(at.rename(columns={'rsid': 'lead_rsid'}), on='lead_rsid', how='left')
es = es.merge(gtf[['gene', 'symbol', 'biotype']], on='gene', how='left')
keep = ['ctx', 'region', 'cell', 'gene', 'symbol', 'biotype', 'num_var', 'lead_rsid', 'start_distance', 'af', 'effect_allele', 'other_allele',
        'allele_status', 'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta', 'qval', 'pval_nominal_threshold', 'is_egene']
save(es[keep], 'egene_stats')
ok('contexts with permutation results', es.ctx.nunique(), 110)
ok('contexts with eGenes', es[es.is_egene].ctx.nunique(), 104)
ok('unique eGenes', es[es.is_egene].gene.nunique(), 3807)

genes = gtf[gtf.gene.isin(es.gene.unique())].copy()
genes['is_egene_anywhere'] = genes.gene.isin(es[es.is_egene].gene)
save(genes, 'genes')
ok('eGenes protein-coding', int(genes[genes.is_egene_anywhere & (genes.biotype == 'protein_coding')].shape[0]), 3517)

# esnp_pairs: pairs of eGenes; independent flag from the clumped official set
m = np.load(B + '/esnp_eg/sets/membership.npz', allow_pickle=True)
eg = np.load(B + '/esnp_eg/sets/membership_eg.npz')['eg']
var = pd.read_csv(B + '/esnp_eg/sets/variants.tsv.gz', sep='\t', usecols=['rsid', 'chr', 'pos', 'MAF'], dtype={'chr': str})
ctx_list = [str(x) for x in m['ctx']]
E = m['e3'] & eg
ind = pd.DataFrame([(var.rsid.iat[i], ctx_list[j]) for i, j in zip(*np.nonzero(E))], columns=['rsid', 'ctx'])
ind['independent'] = True
ok('independent eSNPs (whole brain, unique)', ind.rsid.nunique(), 37412)
pr = pd.read_csv(B + '/esnp_eg/sets/esnp_pairs.tsv.gz', sep='\t')
eg_keys = es.loc[es.is_egene, ['ctx', 'gene']]
pr = pr.merge(eg_keys, on=['ctx', 'gene'], how='inner')
pr = pr.merge(ind, on=['rsid', 'ctx'], how='left')
pr['independent'] = pr.independent.eq(True)
pr = pr.merge(var, on='rsid', how='left').merge(at, on='rsid', how='left')
pr[['region', 'cell']] = pr.ctx.str.split('.', expand=True)
save(pr[['ctx', 'region', 'cell', 'gene', 'biotype', 'rsid', 'chr', 'pos', 'dist', 'MAF', 'af', 'effect_allele', 'other_allele', 'allele_status',
         'pval', 'slope', 'se', 'threshold', 'independent']], 'esnp_pairs')
ok('eSNP-gene pairs of eGenes', len(pr), 3156444)

# contexts
fr = pd.read_csv(B + '/local_inputs/cell_frac_mean.tsv', sep='\t')
n_ind = ind.groupby('ctx').size().rename('independent_esnps')
cx = es.groupby('ctx').agg(genes_tested=('gene', 'size'), egenes=('is_egene', 'sum')).join(n_ind).reset_index()
lnc = es[es.is_egene & (es.biotype != 'protein_coding')].groupby('ctx').size().rename('egenes_lncrna')
cx = cx.merge(lnc, on='ctx', how='left').merge(fr[['ctx', 'region', 'cell', 'n_samples', 'mean_frac', 'median_frac']], on='ctx', how='left')
cx = cx.fillna({'independent_esnps': 0, 'egenes_lncrna': 0})
cx['region'] = pd.Categorical(cx.region, REG); cx['cell'] = pd.Categorical(cx.cell, CELL)
cx = cx.sort_values(['region', 'cell'])
cx['region'] = cx.region.astype(str); cx['cell'] = cx.cell.astype(str)
save(cx[['ctx', 'region', 'cell', 'n_samples', 'mean_frac', 'median_frac', 'genes_tested', 'egenes', 'egenes_lncrna', 'independent_esnps']], 'contexts')

# samples
sa = pd.read_csv(IN + '/database/box-plot/sex-age.csv')
sa['Individual'] = sa.Individual.str.replace('-', '_', regex=False)
n_conf = int((sa.groupby('Individual')[['Sex', 'Age']].nunique() > 1).any(axis=1).sum())
check.append(f"INFO sex-age.csv: {len(sa)} rows, {sa.Individual.nunique()} individuals, individuals with conflicting sex or age: {n_conf}")
sa = sa.drop_duplicates('Individual')
rows = []
for r in REG:
    d = pd.read_csv(B + f'/local_handoff/scbrain2/BayesPrism_InstaPrism_results/estimated_frac/{r}.estimated_frac.txt', sep=r'\s+', quotechar='"')
    d.index.name = 'individual'
    d = d.reset_index(); d.insert(0, 'region', r)
    rows.append(d)
sm = pd.concat(rows, ignore_index=True)
sm['individual'] = sm.individual.astype(str).str.replace('"', '', regex=False)
sm = sm.merge(sa.rename(columns={'Individual': 'individual', 'Sex': 'sex', 'Age': 'age_group'}), on='individual', how='left')
save(sm[['region', 'individual', 'sex', 'age_group'] + CELL], 'samples')
ok('samples', len(sm), 1231)
ok('individuals', sm.individual.nunique(), 134)
check.append(f"INFO samples without sex / age: {int(sm.sex.isna().sum())}")

# fisher
fs = []
for cls, lab in (('nervous', 'neurological'), ('non-nervous', 'non-neurological')):
    x = pd.read_csv(B + f'/fisher_redo/fisher_all_eg_{cls}.tsv', sep='\t')
    x['class'] = lab
    x['disease'] = x.disease.str.replace('@', ' ', regex=False)
    cnt = x[['E.N', 'E.NN', 'NE.N', 'NE.NN']]
    unstable = ((x['E.NN'] < 3) | (x['NE.N'] < 3)) & (x.OR > 1000)
    x['significant'] = (x['q.value'] <= 0.05) & (x.OR > 1) & (x['OR.lower95'] > 1) & (cnt > 0).all(axis=1) & ~unstable
    fs.append(x)
fs = pd.concat(fs, ignore_index=True).drop(columns=['variant'])
fs.columns = [c.replace('.', '_') for c in fs.columns]
save(fs, 'fisher')
ok('Fisher significant (neurological)', int(fs[(fs['class'] == 'neurological') & fs.significant].shape[0]), 857)

# smr with HEIDI and coloc
sr = pd.read_csv(B + '/smr10_fix/smr_all_results.tsv.gz', sep='\t')
co = pd.read_csv(B + '/coloc10_fix/coloc_all.tsv', sep='\t').rename(columns={'ctx': 'Region-Cell', 'gwas': 'GWAS', 'gene': 'probeID'})
sr = sr.merge(co[['Region-Cell', 'GWAS', 'probeID', 'nsnps', 'PP.H3.abf', 'PP.H4.abf']], on=['Region-Cell', 'GWAS', 'probeID'], how='left')
dz = pd.read_csv(B + '/integrate10_fix/gene_evidence.tsv', sep='\t', usecols=['disease', 'GWAS']).drop_duplicates()
sr = sr.merge(dz, on='GWAS', how='left')
sr['ctx'] = sr['Region-Cell'].str.replace('-', '.', regex=False)
sr = sr.drop(columns=['Region-Cell'])
sr['smr_heidi_coloc'] = (sr.q_SMR < 0.05) & (sr.p_HEIDI >= 0.01) & (sr['PP.H4.abf'] >= 0.8)
sr.columns = [c.replace('.', '_') for c in sr.columns]
save(sr, 'smr')
p3 = sr[sr.smr_heidi_coloc][['disease', 'probeID']].drop_duplicates()
ok('disease-gene pairs with SMR q<0.05', int(sr[sr.q_SMR < 0.05][['disease', 'probeID']].drop_duplicates().shape[0]), 253)
ok('disease-gene pairs with SMR + HEIDI + PP.H4>=0.8', len(p3), 32)

ge = pd.read_csv(B + '/integrate10_fix/gene_evidence.tsv', sep='\t')
save(ge, 'gene_evidence')
sg = pd.read_csv(B + '/scgwas10/scgwas_sig_modules_new.tsv', sep='\t')
save(sg, 'scgwas_modules')
ok('scGWAS significant modules', len(sg), 118)

open(OUT + '/CHECK.txt', 'w').write('\n'.join(check) + '\n')
print('\n'.join(check))
print('S22_DONE')
