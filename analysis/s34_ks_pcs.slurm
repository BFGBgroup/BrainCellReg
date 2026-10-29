#!/bin/bash
# KS of the 20 Enformer PCs by group (s34_ks_pcs.py), PI #84.
#SBATCH -J scb_kspc
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s34_kspc_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s34_kspc_%j.err
source /etc/profile
cd /public/home/yucong/scbrain/scripts
/public/home/yucong/.conda/envs/scvi/bin/python s34_ks_pcs.py || echo KS_PCS_FAILED
