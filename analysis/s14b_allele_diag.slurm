#!/bin/bash
#SBATCH -J scb_adiag
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH -t 00:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s14b_adiag_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s14b_adiag_%j.err
# 单核复现的方向一致率低于 50%，查两边的等位基因标签（对 1000G 欧洲参考）
OUT=/public/home/yucong/scbrain/sn_replication
G=/public/home/yucong/scbrain/local_handoff/scbrain1/g1000_eur/g1000_eur
/public/home/yucong/scbrain/bin/plink --bfile $G --freq --out $OUT/g1000 > /dev/null || { echo PLINK_FAILED; exit 1; }
/public/home/yucong/.conda/envs/scvi/bin/python /public/home/yucong/scbrain/scripts/s14b_allele_diag.py
