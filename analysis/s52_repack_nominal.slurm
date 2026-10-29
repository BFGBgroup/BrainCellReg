#!/bin/bash
#SBATCH -J s52_repack
#SBATCH -p debug
#SBATCH -n 16
#SBATCH --mem=80G
#SBATCH -t 00:59:00
#SBATCH -o /public/home/yucong/scbrain/logs/s52_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s52_%j.err
# debug 每人只能提交有限个作业，所以一个作业里并行 8 个 worker；57 分钟到点自停，靠 _log 续跑，重复提交直到 2,420 个都 ok
source /etc/profile
module load apps/anaconda3/2024.02
for k in 0 1 2 3 4 5 6 7; do
  timeout 57m python /public/home/yucong/scbrain/scripts/s52_repack_nominal.py $k 8 &
done
wait
# 变异索引（s53，几分钟）；debug 每人只能有 1 个作业，所以搭在这里，已建好就跳过
[ -f /public/home/yucong/scbrain/db/tables/variants.parquet ] || python /public/home/yucong/scbrain/scripts/s53_variants.py
cat /public/home/yucong/scbrain/nominal_web/_log/task*.tsv | awk -F'\t' '$7=="ok"' | wc -l
