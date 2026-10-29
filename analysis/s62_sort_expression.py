"""网站提速（韩教授 10-20 #146「用起来有点卡…是检索的问题」）：把网站查询用到的表按查询键重排、写小 row group，
让 DuckDB 按 row group 的 min/max 跳过无关部分（原表无序：查一个基因要扫 expression 全部 4,143 个 row group，3.3 s）。
同名表内容不变（逐表核对行数 + 全行哈希，不等就报错，不写 DONE）；另出 3 份按 rsid 排序的副本给变异页用。
in   db/tables/{expression,smr,expr_summary}.parquet、db/tables_allele_v3/{egene_stats,esnp_pairs,geno_box}.parquet（与线上 md5 相同）
out  db/tables_sorted/<同名>.parquet + esnp_pairs_by_rsid / geno_box_by_rsid / egene_stats_by_lead.parquet"""
import os, time
import duckdb

B = '/public/home/yucong/scbrain/db'
OUT = f'{B}/tables_sorted'
os.makedirs(OUT, exist_ok=True)
JOBS = [  # (output name, source, ORDER BY)
    ('expression', f'{B}/tables/expression.parquet', 'gene, ctx, sample'),
    ('egene_stats', f'{B}/tables_allele_v3/egene_stats.parquet', 'gene, ctx'),
    ('esnp_pairs', f'{B}/tables_allele_v3/esnp_pairs.parquet', 'gene, ctx, pval'),
    ('geno_box', f'{B}/tables_allele_v3/geno_box.parquet', 'gene, rsid, ctx'),
    ('expr_summary', f'{B}/tables/expr_summary.parquet', 'gene, ctx'),
    ('smr', f'{B}/tables/smr.parquet', 'probeID, disease, ctx'),
    ('esnp_pairs_by_rsid', f'{B}/tables_allele_v3/esnp_pairs.parquet', 'rsid, gene, ctx'),
    ('geno_box_by_rsid', f'{B}/tables_allele_v3/geno_box.parquet', 'rsid, gene, ctx'),
    ('egene_stats_by_lead', f'{B}/tables_allele_v3/egene_stats.parquet', 'lead_rsid, gene, ctx'),
]
con = duckdb.connect()
con.execute("SET threads TO 8"); con.execute("SET memory_limit = '50GB'"); con.execute(f"SET temp_directory = '{OUT}/tmp'")
print('duckdb', duckdb.__version__, flush=True)
for name, src, order in JOBS:
    dst = f'{OUT}/{name}.parquet'
    t = time.time()
    con.execute(f"COPY (SELECT * FROM read_parquet('{src}') ORDER BY {order}) TO '{dst}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 20000)")
    cols = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{src}')").fetchall()]
    h = f"SELECT count(*), sum(hash({', '.join(chr(34) + c + chr(34) for c in cols)})::HUGEINT) FROM read_parquet('{{}}')"
    a = con.execute(h.format(src)).fetchone(); b = con.execute(h.format(dst)).fetchone()
    same_cols = cols == [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{dst}')").fetchall()]
    rg = con.execute(f"SELECT count(DISTINCT row_group_id) FROM parquet_metadata('{dst}')").fetchone()[0]
    print(f'{name:20s} {time.time() - t:5.0f} s  {os.path.getsize(dst) / 1e6:6.0f} MB  rows {a[0]:,}  row groups {rg}  '
          f'same rows+hash {a == b}  same columns {same_cols}', flush=True)
    assert a == b and same_cols, f'CONTENT DIFFERS: {name}'
t = time.time()
n = con.execute(f"SELECT count(*) FROM read_parquet('{OUT}/expression.parquet') WHERE gene = 'ENSG00000136235' AND ctx = 'FCTX.ExN'").fetchone()[0]
print(f'check: GPNMB FCTX.ExN expression rows {n} in {time.time() - t:.2f} s')
print('S62_DONE')
