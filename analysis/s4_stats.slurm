#!/bin/bash
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --mem=100G
#SBATCH -t 59:00
#SBATCH -J scb_stats
#SBATCH -o /public/home/yucong/scbrain/logs/%x_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/%x_%j.err

# 第 5–10、12 步 + MAF 分箱的统计表
source /etc/profile
export LC_ALL=C
export NPROC=32
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /public/home/yucong/scbrain
export SCB_OUT=${SCB_OUT:-/public/home/yucong/scbrain/esnp_func}
echo "SCB_OUT=$SCB_OUT ESNP_DIR=${ESNP_DIR:-default}"
date
[ -f $SCB_OUT/sets/membership_eg.npz ] || /public/home/yucong/.conda/envs/scvi/bin/python scripts/s1b_egene_flag.py || exit 1
date
/public/home/yucong/.conda/envs/scvi/bin/python scripts/s4_stats.py
echo "exit=$?"
date
ls -la $SCB_OUT/results
