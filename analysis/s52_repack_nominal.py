"""全量 cis nominal 打包成网站用的格式（s51 试打包定的规格）。
输入 paper2 .../run/new/<ctx>.cis_qtl_pairs.<chr>.parquet（2,420 个，只读）
输出 OUT/<chr>/<ctx>.parquet：pos(int32, bim 第 4 列，同 s26 的 pos38) phenotype_id variant_id af slope slope_se mlog10p(float32)
按 (pos, phenotype_id) 排序，zstd 9，行组 25 万行；每个文件写完读回逐值核对，结果记 OUT/_log/task<k>.tsv。
用法：python s52_repack_nominal.py <task_id> <n_tasks>；已存在且核对过的跳过（可重跑）。"""
import glob, os, sys, time
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

IN = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data/'
D = IN + 'run/new/'
OUT = '/public/home/yucong/scbrain/nominal_web/'
k, n = int(sys.argv[1]), int(sys.argv[2])
fs = sorted(glob.glob(D + '*.cis_qtl_pairs.chr*.parquet'))
assert len(fs) == 2420, len(fs)
mine = fs[k::n]
os.makedirs(OUT + '_log', exist_ok=True)
logf = OUT + f'_log/task{k}.tsv'
done = set()
if os.path.exists(logf):
    done = {ln.split('\t')[0] for ln in open(logf) if ln.rstrip().endswith('\tok')}

bim = pd.read_csv(IN + 'genotype/out.bim', sep='\t', header=None, usecols=[1, 3], names=['rsid', 'pos'],
                  dtype={'rsid': str, 'pos': np.int64}).drop_duplicates('rsid')
pos_of = pd.Series(bim.pos.values, index=bim.rsid.values)

for src in mine:
    name = os.path.basename(src)
    if name in done:
        continue
    t0 = time.time()
    ctx, chrom = name.split('.cis_qtl_pairs.')[0], name.split('.cis_qtl_pairs.')[1].replace('.parquet', '')
    t = pq.read_table(src, columns=['phenotype_id', 'variant_id', 'af', 'pval_nominal', 'slope', 'slope_se']).to_pandas()
    pos = pos_of.reindex(t.variant_id.values).values
    n_nopos = int(np.isnan(pos).sum())
    t['pos'] = np.where(np.isnan(pos), -1, pos).astype(np.int32)
    t['mlog10p'] = (-np.log10(np.maximum(t.pval_nominal.values, 1e-300))).astype(np.float32)
    t = t.sort_values(['pos', 'phenotype_id'], kind='stable').reset_index(drop=True)
    nt = pa.table({
        'pos': pa.array(t.pos.values),
        'phenotype_id': pa.array(t.phenotype_id.values).dictionary_encode(),
        'variant_id': pa.array(t.variant_id.values),
        'af': pa.array(t.af.values.astype(np.float32)),
        'slope': pa.array(t.slope.values.astype(np.float32)),
        'slope_se': pa.array(t.slope_se.values.astype(np.float32)),
        'mlog10p': pa.array(t.mlog10p.values),
    })
    os.makedirs(OUT + chrom, exist_ok=True)
    out = OUT + f'{chrom}/{ctx}.parquet'
    pq.write_table(nt, out + '.tmp', compression='zstd', compression_level=9, row_group_size=250_000)
    back = pq.read_table(out + '.tmp').to_pandas()
    ok = (len(back) == len(t) == pq.ParquetFile(src).metadata.num_rows
          and (back.variant_id.values == t.variant_id.values).all()
          and (back.phenotype_id.astype(str).values == t.phenotype_id.values).all()
          and np.allclose(back.slope, t.slope, equal_nan=True)
          and np.allclose(back.slope_se, t.slope_se, equal_nan=True)
          and np.allclose(10 ** (-back.mlog10p.astype(float)), t.pval_nominal, rtol=1e-5, atol=0, equal_nan=True))
    if ok:
        os.replace(out + '.tmp', out)
    with open(logf, 'a') as fh:
        fh.write(f'{name}\t{len(t)}\t{n_nopos}\t{os.path.getsize(src)}\t{os.path.getsize(out if ok else out + ".tmp")}'
                 f'\t{time.time() - t0:.0f}\t{"ok" if ok else "MISMATCH"}\n')
    print(name, len(t), 'ok' if ok else 'MISMATCH', flush=True)
