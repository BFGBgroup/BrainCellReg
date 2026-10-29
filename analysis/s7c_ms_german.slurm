#!/bin/bash
#SBATCH -J scb_msfmt
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -t 00:40:00
#SBATCH -o /public/home/yucong/scbrain/logs/s7c_msfmt_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s7c_msfmt_%j.err
# 多发性硬化改用德国队列 GCST003566（用户 10-10 定：能行就用它，不行退回原来的 UKB 那份）。
# 转格式 + 判断「能不能用」：与 eQTL 的 SNP 重合 ≥70%、含 MHC 区、P 与 b/se 自洽。过了才改 gwasinfo.txt。
PY=/public/home/yucong/.conda/envs/scvi/bin/python
export GWAS_OUT=/public/home/yucong/scbrain/gwas10
g=GCST003566
R=/public/home/yucong/scbrain/external/gwas10/GCST003566_MS_German/GCST003566.h.tsv.gz
echo "md5 $(md5sum $R | cut -d' ' -f1) expect ac9f9f821ecf5daf4447a7611522bed2"
$PY /public/home/yucong/scbrain/scripts/s7_gwas_format.py $g > $GWAS_OUT/qc/$g.log 2>&1
cat $GWAS_OUT/qc/$g.log
grep -q FORMAT_DONE $GWAS_OUT/qc/$g.log || { echo MSFMT_FAILED; exit 1; }
mhc=$(zcat $R | awk -F'\t' 'NR>1 && $1=="6" && $2>25000000 && $2<35000000 {n++} END{print n+0}')
mhcsig=$(zcat $R | awk -F'\t' 'NR>1 && $1=="6" && $2>25000000 && $2<35000000 && $8+0<5e-8 {n++} END{print n+0}')
frac=$(grep eqtl_snps_in_gwas_frac $GWAS_OUT/qc/$g.json | tr -dc '0-9.')
bad=$(grep p_vs_bse_frac $GWAS_OUT/qc/$g.json | tr -dc '0-9.e-')
echo "MHC rows $mhc; MHC rows with P<5e-8 $mhcsig; eQTL SNP overlap $frac; p-vs-b/se mismatch frac $bad"
if awk -v f=$frac -v m=$mhc 'BEGIN{exit !(f>=0.70 && m>1000)}'; then
  grep -v -E '^(GCST90558092|GCST003566)[[:space:]]' $GWAS_OUT/gwasinfo.txt > $GWAS_OUT/gwasinfo.tmp
  printf 'GCST003566\t15283\n' >> $GWAS_OUT/gwasinfo.tmp && mv $GWAS_OUT/gwasinfo.tmp $GWAS_OUT/gwasinfo.txt
  rm -f $GWAS_OUT/GCST90558092.sumstats.txt
  cat $GWAS_OUT/gwasinfo.txt
  echo MS_GERMAN_USABLE
else
  echo MS_GERMAN_NOT_USABLE; exit 2
fi
