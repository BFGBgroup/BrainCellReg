#!/bin/bash
#SBATCH -J scb_gwasfmt
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 01:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s7_gwasfmt_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s7_gwasfmt_%j.err
# 7 份新 GWAS 转 SMR / coloc / scGWAS 格式；没换的 3 份直接链到贺奕捷转好的文件
PY=/public/home/yucong/.conda/envs/scvi/bin/python
S=/public/home/yucong/scbrain/scripts
IN=/public/home/yucong/scbrain_in
export GWAS_OUT=/public/home/yucong/scbrain/gwas10
mkdir -p $GWAS_OUT/qc
for g in GCST90432156 GCST90558092 GCST90726344 GCST009325 GCST90018844 GCST90128471 GCST011179; do
  $PY $S/s7_gwas_format.py $g > $GWAS_OUT/qc/$g.log 2>&1 &
done
wait
for g in GCST013196 GCST90014122 GCST90477562; do
  ln -sf $IN/paper2/10disease-gwas/$g-smr.tsv $GWAS_OUT/$g-smr.tsv
  ln -sf $IN/paper2/10disease-gwas/$g-coloc.tsv $GWAS_OUT/$g-coloc.tsv
  ln -sf $IN/scGWAS_inputs/run/GWAS/$g.sumstats.txt $GWAS_OUT/$g.sumstats.txt
done
# MAGMA 用的总样本量（病例 + 对照）；没换的 3 份沿用他的 gwasinfo.txt
printf 'GCST90432156\t437736\nGCST90558092\t360004\nGCST90726344\t514455\nGCST009325\t482730\nGCST90018844\t521177\nGCST90128471\t130644\nGCST011179\t357806\n' > $GWAS_OUT/gwasinfo.txt
grep -E '^(GCST013196|GCST90014122|GCST90477562)[[:space:]]' $IN/scGWAS_inputs/run/GWAS/gwasinfo.txt >> $GWAS_OUT/gwasinfo.txt
cat $GWAS_OUT/gwasinfo.txt
tail -n 3 $GWAS_OUT/qc/*.log
ls -la $GWAS_OUT | cut -c25-
n=$(grep -l FORMAT_DONE $GWAS_OUT/qc/*.log | wc -l)
echo "formatted $n / 7"
[ "$n" = 7 ] && echo GWASFMT_DONE || echo GWASFMT_INCOMPLETE
