#!/bin/bash
# Candidate gene-set sizes for the Fig 6 GO panel (count_gwas_egenes.py), PI #85.
#SBATCH -J scb_gocount
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH -o /public/home/yucong/scbrain/logs/gocount_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/gocount_%j.err
source /etc/profile
cd /public/home/yucong/scbrain/scripts
/public/home/yucong/.conda/envs/scvi/bin/python -u count_gwas_egenes.py || echo GOCOUNT_FAILED
