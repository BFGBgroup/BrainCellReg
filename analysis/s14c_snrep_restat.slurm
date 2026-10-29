#!/bin/bash
#SBATCH -J scb_snrep2
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s14c_snrep2_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s14c_snrep2_%j.err
# 复现检验重算统计：等位基因朝向按 1000G 欧洲频率逐位点判定。第一遍（按 bim 字母）的结果留在 run1_labels_as_given/
OUT=/public/home/yucong/scbrain/sn_replication
mkdir -p $OUT/run1_labels_as_given
for f in rep_by_context.tsv rep_matrix.tsv criteria.txt align_qc.tsv pi1_8x8_main_FCTX_TCTX_WHMT.tsv pi1_8x8_other_7_regions.tsv pi1_8x8_all_10_regions.tsv; do
  [ -f $OUT/$f ] && [ ! -f $OUT/run1_labels_as_given/$f ] && cp -p $OUT/$f $OUT/run1_labels_as_given/$f
done
/public/home/yucong/.conda/envs/scvi/bin/python /public/home/yucong/scbrain/scripts/s14_snrep.py stats || { echo STATS_FAILED; exit 1; }
echo "RESTAT_DONE $(date '+%F %T')"
