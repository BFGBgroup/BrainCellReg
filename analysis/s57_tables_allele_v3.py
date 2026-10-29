"""Website tables with the allele_table_v3 labels (PI 10-19 #135, option A). Only the letters change: the counted allele of the
dosages is the same, v3 just names it for variants v2 left unresolved; n_effect in geno_box (copies of the counted allele) stays.
  in   db/tables/{variants,egene_stats,esnp_pairs}.parquet, db/geno_box_all.parquet, allele_fix/allele_table_v{2,3}.tsv.gz
  out  db/tables_allele_v3/{variants,egene_stats,esnp_pairs,geno_box}.parquet (old tables untouched)
       variants gets in_1000g (v2 status != not_in_1000G): the website's "not clumped" label, which v3 statuses no longer carry.
  check  rows unchanged; letters of rows already resolved in v2 unchanged; counts of resolved rows before / after."""
import os, sys
import duckdb
import pandas as pd

B = '/public/home/yucong/scbrain'
T = B + '/db/tables'
OUT = B + '/db/tables_allele_v3'
os.makedirs(OUT, exist_ok=True)
RES = ['ok', 'ok_ld', 'ok_strand', 'ok_alfa', 'ok_alfa_strand']
con = duckdb.connect()

v3 = pd.read_csv(B + '/allele_fix/allele_table_v3.tsv.gz', sep='\t', usecols=['snp', 'eff', 'oth', 'status']).drop_duplicates('snp')
v3.loc[~v3.status.isin(RES), ['eff', 'oth']] = None
v2 = pd.read_csv(B + '/allele_fix/allele_table_v2.tsv.gz', sep='\t', usecols=['snp', 'status']).drop_duplicates('snp')
v3 = v3.merge(v2.rename(columns={'status': 'status_v2'}), on='snp', how='left')
v3['in_1000g'] = v3.status_v2 != 'not_in_1000G'
v3 = v3.rename(columns={'eff': 'effect_allele', 'oth': 'other_allele', 'status': 'allele_status'})
bad = []

for name, src, key in (('variants', T + '/variants.parquet', 'rsid'), ('egene_stats', T + '/egene_stats.parquet', 'lead_rsid'),
                       ('esnp_pairs', T + '/esnp_pairs.parquet', 'rsid'), ('geno_box', B + '/db/geno_box_all.parquet', 'rsid')):
    d = con.execute(f"SELECT * FROM read_parquet('{src}')").df()
    cols = list(d.columns)
    old = d[['effect_allele', 'other_allele', 'allele_status']].copy()
    lab = v3[['snp', 'effect_allele', 'other_allele', 'allele_status'] + (['in_1000g'] if name == 'variants' else [])]
    m = d.drop(columns=['effect_allele', 'other_allele', 'allele_status']).merge(
        lab.rename(columns={'snp': key}), on=key, how='left')
    if name == 'variants':
        cols = cols + ['in_1000g']
    m = m[cols]
    n_before, n_after = int(old.effect_allele.notna().sum()), int(m.effect_allele.notna().sum())
    was = old.allele_status.isin(['ok', 'ok_ld']).values
    same = ((old.effect_allele.values[was] == m.effect_allele.values[was]) & (old.other_allele.values[was] == m.other_allele.values[was])).mean()
    print(f'{name}: {len(d):,} rows -> {len(m):,}; with effect allele {n_before:,} ({n_before / len(d):.1%}) -> {n_after:,} ({n_after / len(m):.1%}); '
          f'v2-resolved rows with same letters {same:.6f}', flush=True)
    print(m.allele_status.value_counts(dropna=False).to_string(), flush=True)
    if len(m) != len(d) or same < 1:
        bad.append(name)
    con.register('m', m)
    order = ' ORDER BY rsid' if name == 'variants' else ''
    rg = ', ROW_GROUP_SIZE 100000' if name == 'variants' else ''
    con.execute(f"COPY (SELECT * FROM m{order}) TO '{OUT}/{name}.parquet' (FORMAT PARQUET, COMPRESSION ZSTD{rg})")
    con.unregister('m')

if bad:
    sys.exit('CHECK_FAILED ' + ' '.join(bad))
print('S57_DONE')
