#!/bin/bash
# Per-combination features for the Fig 4g atlas (s32_ctx_features.py).
#SBATCH -J scb_ctxfeat
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s32_ctxfeat_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s32_ctxfeat_%j.err
source /etc/profile
module load apps/bedtools/2.30.0
cd /public/home/yucong/scbrain/scripts
/public/home/yucong/.conda/envs/scvi/bin/python s32_ctx_features.py || echo CTX_FEATURES_FAILED
