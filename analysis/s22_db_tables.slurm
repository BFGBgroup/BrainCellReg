#!/bin/bash
#SBATCH -J scb_dbtab
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s22_dbtab_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s22_dbtab_%j.err
/public/home/yucong/.conda/envs/scvi/bin/python /public/home/yucong/scbrain/scripts/s22_db_tables.py
