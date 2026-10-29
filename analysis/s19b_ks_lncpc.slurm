#!/bin/bash
# KS / ECDF for the lncRNA and protein-coding eSNP subsets (Fig 5, PI 10-14 #33 / #36)
#SBATCH -J scb_ks2
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 00:40:00
#SBATCH -o /public/home/yucong/scbrain/logs/s19b_ks_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s19b_ks_%j.err
cd /public/home/yucong/scbrain
export SCB_OUT=/public/home/yucong/scbrain/esnp_eg NPROC=8
for s in lnc pc; do ESET=$s /public/home/yucong/.conda/envs/scvi/bin/python scripts/s19_ks_ctx.py || exit 1; done
echo "S19B_DONE"
