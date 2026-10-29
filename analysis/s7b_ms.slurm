#!/bin/bash
#SBATCH -J scb_msfmt
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -t 00:40:00
#SBATCH -o /public/home/yucong/scbrain/logs/s7b_msfmt_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s7b_msfmt_%j.err
# 多发性硬化单独重转：用作者原始文件，按位置到 1000G 欧洲参考里取 rsID
PY=/public/home/yucong/.conda/envs/scvi/bin/python
export GWAS_OUT=/public/home/yucong/scbrain/gwas10
R=/public/home/yucong/scbrain/external/gwas10/GCST90558092_MS/GCST90558092.tsv
echo "md5 $(md5sum $R | cut -d' ' -f1) expect 2a00fe3fabf21631bf62a7e62699ba4d"
$PY /public/home/yucong/scbrain/scripts/s7_gwas_format.py GCST90558092 > $GWAS_OUT/qc/GCST90558092.log 2>&1
cat $GWAS_OUT/qc/GCST90558092.log
grep -q FORMAT_DONE $GWAS_OUT/qc/GCST90558092.log && echo MSFMT_DONE || echo MSFMT_FAILED
