#!/bin/bash
# 单核复现检验敏感性：等位基因朝向改用 allele_table_v2（频率 + 连锁判定），判据不变；结果写 sn_replication/run3_allele_v2/，正式结果不动
#SBATCH -J scb_snrep3
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 00:45:00
#SBATCH -o /public/home/yucong/scbrain/logs/s14d_snrep3_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s14d_snrep3_%j.err
export ALLELE_TABLE=/public/home/yucong/scbrain/allele_fix/allele_table_v2.tsv.gz
export SNREP_RES=/public/home/yucong/scbrain/sn_replication/run3_allele_v2
/public/home/yucong/.conda/envs/scvi/bin/python /public/home/yucong/scbrain/scripts/s14_snrep.py stats || { echo STATS_FAILED; exit 1; }
echo "SNREP3_DONE $(date '+%F %T')"
