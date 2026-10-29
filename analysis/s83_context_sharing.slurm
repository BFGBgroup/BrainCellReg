#!/bin/bash
#SBATCH -J scb_sharing
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH -o /public/home/yucong/scbrain/logs/fig3_prep_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/fig3_prep_%j.err
mkdir -p /public/home/yucong/scbrain/figures/fig3
/public/home/yucong/.conda/envs/scvi/bin/python /public/home/yucong/scbrain/scripts/s83_context_sharing.py
