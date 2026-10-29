#!/bin/bash
#SBATCH -J scb_ks
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 01:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s19_ks_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s19_ks_%j.err
cd /public/home/yucong/scbrain
export SCB_OUT=/public/home/yucong/scbrain/esnp_eg NPROC=8
/public/home/yucong/.conda/envs/scvi/bin/python scripts/s19_ks_ctx.py
