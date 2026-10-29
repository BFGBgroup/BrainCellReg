"""Locus plots for every gene (PI 10-19「这里的疾病无法选，出图后选项就消失了」): GWAS P of the 10 diseases and an LD panel for
every variant of the site, not only the 4 demo genes.
1. gwas_web/<disease>.parquet: rsid, chr (chr1..), pos (bim, = s26 pos38), mlog10p (float32) for the variants in out.bim,
   from gwas10/<acc>-smr.tsv (the files SMR used), sorted by chr, pos.
2. g1000_web/chr<N>.bed/bim/fam: 1000 Genomes EUR (503 people, the panel of s26's demo LD), cut to the variants of out.bim
   and split by chromosome, so the site computes r2 with plink on demand for any index variant.
Checks go to db/CHECK_web_gwas_ld.txt."""
import os, subprocess
import numpy as np
import pandas as pd
import duckdb

B = '/public/home/yucong/scbrain'
IN = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data'
OUT_G = B + '/db/gwas_web'
OUT_L = B + '/db/g1000_web'
os.makedirs(OUT_G, exist_ok=True); os.makedirs(OUT_L, exist_ok=True)
DZ = [('PD', 'GCST009325'), ('AD', 'GCST013196'), ('MS', 'GCST003566'), ('MG', 'GCST90432156'), ('SCZ', 'GCST90128471'),
      ('MDD', 'GCST90726344'), ('Narcolepsy', 'GCST90477562'), ('LacunarStroke', 'GCST90014122'), ('CUD', 'GCST011179'),
      ('FS', 'GCST90018844')]
check = []
bim = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[0, 1, 3], names=['chr', 'rsid', 'pos'],
                  dtype={'chr': str}).drop_duplicates('rsid')
con = duckdb.connect()
for dz, acc in DZ:
    f = f'{B}/gwas10/{acc}-smr.tsv'
    x = pd.read_csv(f, sep='\t', usecols=['SNP', 'p'], dtype={'SNP': str})
    n0 = len(x)
    x = x.rename(columns={'SNP': 'rsid'}).drop_duplicates('rsid').merge(bim, on='rsid', how='inner')
    p = pd.to_numeric(x.p, errors='coerce').values
    x['mlog10p'] = (-np.log10(np.clip(p, 1e-300, 1))).astype(np.float32)
    x = x[np.isfinite(x.mlog10p)]
    x['pos'] = x.pos.astype(np.int32)
    con.register('x', x[['rsid', 'chr', 'pos', 'mlog10p']])
    con.execute(f"COPY (SELECT * FROM x ORDER BY chr, pos) TO '{OUT_G}/{dz}.parquet' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)")
    con.unregister('x')
    check.append(f'GWAS {dz} {acc}: {n0:,} rows in file -> {len(x):,} on out.bim variants; min P {10 ** -float(x.mlog10p.max()):.3g}')
    print(check[-1], flush=True)

# 1000G EUR cut to our variants, per chromosome
bf = B + '/local_handoff/scbrain1/g1000_eur/g1000_eur'
snps = OUT_L + '/variants.txt'
bim.rsid.to_csv(snps, index=False, header=False)
for c in [str(i) for i in range(1, 23)]:
    r = subprocess.run([B + '/bin/plink', '--bfile', bf, '--extract', snps, '--chr', c, '--make-bed', '--out', f'{OUT_L}/chr{c}',
                        '--threads', '4', '--memory', '8000'], capture_output=True, text=True)
    nb = sum(1 for _ in open(f'{OUT_L}/chr{c}.bim')) if os.path.exists(f'{OUT_L}/chr{c}.bim') else -1
    check.append(f'1000G chr{c}: {nb:,} variants (plink exit {r.returncode})')
    print(check[-1], flush=True)
for f in os.listdir(OUT_L):
    if f.endswith(('.log', '.nosex')):
        os.remove(os.path.join(OUT_L, f))
open(B + '/db/CHECK_web_gwas_ld.txt', 'w').write('\n'.join(check) + '\n')
