"""Candidate gene sets for the Fig 6 GO panel (PI #85), counted before any enrichment is run:
  A  SMR genes: q_SMR < 0.05 in >= 1 context (integrate10_fix/gene_evidence.tsv; as s20_enrich.R)
  B  GWAS-locus eGenes: eGenes (q < 0.05 in >= 1 context) with >= 1 eSNP at genome-wide significance (P < 5e-8) in the GWAS
Prints genes per disease, and how many are in the MHC (chr6 25-35 Mb)."""
import duckdb, pandas as pd
B = '/public/home/yucong/scbrain'
T = B + '/db/tables'
NAME = {'GCST90432156': 'MG', 'GCST003566': 'MS', 'GCST90726344': 'DEP', 'GCST009325': 'PD', 'GCST90018844': 'FS',
        'GCST90128471': 'SCZ', 'GCST011179': 'CD', 'GCST013196': 'AD', 'GCST90014122': 'LS', 'GCST90477562': 'NAR'}
c = duckdb.connect()
c.execute("SET threads TO 4")
ge = pd.read_csv(B + '/integrate10_fix/gene_evidence.tsv', sep='\t')
print('A  SMR genes per disease:', ge.groupby('disease').size().to_dict(), '| union', ge.iloc[:, 0].nunique() if len(ge) else 0)
c.execute(f"""CREATE TABLE eg AS SELECT DISTINCT p.rsid, p.gene FROM '{T}/esnp_pairs.parquet' p
              JOIN (SELECT DISTINCT gene FROM '{T}/egene_stats.parquet' WHERE qval < 0.05) e USING (gene)""")
print('eGene eSNP pairs', c.execute("SELECT count(*), count(DISTINCT gene) FROM eg").fetchone())
genes = c.execute(f"SELECT gene, chr, start FROM '{T}/genes.parquet'").df()
rows = []
for g, dz in NAME.items():
    n = c.execute(f"""SELECT DISTINCT eg.gene FROM eg JOIN read_csv_auto('{B}/gwas10/{g}.sumstats.txt', delim='\t', header=true) s
                      ON lower(s.SNP) = eg.rsid WHERE s.P < 5e-8""").df()
    m = n.merge(genes, on='gene', how='left')
    mhc = int(((m.chr.astype(str) == '6') & (m.start > 25e6) & (m.start < 35e6)).sum())
    rows.append(dict(disease=dz, B_genes=len(n), B_mhc=mhc))
print(pd.DataFrame(rows).to_string())
