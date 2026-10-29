#!/bin/bash
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem=60G
#SBATCH -t 1:00:00
#SBATCH -J scb_fisher
#SBATCH -o /public/home/yucong/scbrain/logs/%x_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/%x_%j.err

# 疾病风险位点富集重跑（eSNP 改为只留 eGene 的），同时对回原结果
source /etc/profile
module load apps/R/4.5.0
export LC_ALL=C NPROC=22
cd /public/home/yucong/scbrain
date
/public/home/yucong/.conda/envs/scvi/bin/python scripts/s6_fisher_counts.py | tail -4 || exit 1
date
Rscript scripts/s6_fisher_test.R || exit 1
date
ls -la fisher_redo | head -30
