"""只读元数据：全量 nominal parquet 的列、行数、行组、压缩（不读数据，登录节点可跑）。"""
import glob, os
import pyarrow.parquet as pq

D = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data/run/new/'
fs = sorted(glob.glob(D + '*.cis_qtl_pairs.chr*.parquet'))
print('files', len(fs), 'bytes', sum(os.path.getsize(f) for f in fs))
m = pq.ParquetFile(fs[0]).metadata
print(fs[0], m.num_rows, 'rows', m.num_row_groups, 'rg')
print(pq.ParquetFile(fs[0]).schema_arrow)
rg = m.row_group(0)
for i in range(rg.num_columns):
    c = rg.column(i)
    print(c.path_in_schema, c.compression, c.total_compressed_size, c.total_uncompressed_size)
n = 0
for f in fs:
    n += pq.ParquetFile(f).metadata.num_rows
print('total rows', n)
print(pq.read_table(fs[0]).slice(0, 5).to_pandas().to_string())
