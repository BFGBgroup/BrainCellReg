"""全部检验过的变异的索引：rs 号 -> 染色体、位置（bim，同 s26 pos38）、等位基因（allele_table_v2，同 s22）。
网站据此找到全量 nominal 的文件（nominal_web/<chr>/）和行（pos）。输出 db/tables/variants.parquet，按 rsid 排序。"""
import duckdb, pandas as pd

B = '/public/home/yucong/scbrain'
IN = '/public/home/yucong/scbrain_in/paper2'
bim = pd.read_csv(IN + '/single-cell-eqtl/data/genotype/out.bim', sep='\t', header=None, usecols=[0, 1, 3, 4, 5],
                  names=['chr', 'rsid', 'pos', 'a1', 'a2'], dtype={'chr': str})
n0 = len(bim)
dup = bim.rsid.duplicated(keep=False).sum()
bim = bim.drop_duplicates('rsid')
at = pd.read_csv(B + '/allele_fix/allele_table_v2.tsv.gz', sep='\t', usecols=['snp', 'eff', 'oth', 'status'])
at = at.drop_duplicates('snp').rename(columns={'snp': 'rsid', 'eff': 'effect_allele', 'oth': 'other_allele', 'status': 'allele_status'})
at.loc[~at.allele_status.isin(['ok', 'ok_ld']), ['effect_allele', 'other_allele']] = None
v = bim.merge(at, on='rsid', how='left').sort_values('rsid').reset_index(drop=True)
v['pos'] = v.pos.astype('int32')
con = duckdb.connect()
con.register('v', v)
con.execute(f"COPY (SELECT * FROM v ORDER BY rsid) TO '{B}/db/tables/variants.parquet' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)")
print('bim rows', n0, 'rows with duplicated rsid', dup, 'variants', len(v))
print('with allele table', int(v.allele_status.notna().sum()), 'ok/ok_ld', int(v.effect_allele.notna().sum()))
print(v.chr.value_counts().sort_index().to_string())
