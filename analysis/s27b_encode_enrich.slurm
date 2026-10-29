#!/bin/bash
# ENCODE prefrontal-cortex histone / DNase enrichment of eSNPs + MAF by group (s27b_encode_enrich.py).
#SBATCH -J scb_encode
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s27b_encode_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s27b_encode_%j.err
source /etc/profile
module load apps/bedtools/2.30.0
cd /public/home/yucong/scbrain/scripts
/public/home/yucong/.conda/envs/scvi/bin/python s27b_encode_enrich.py || echo ENCODE_ENRICH_FAILED
