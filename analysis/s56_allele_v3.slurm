#!/bin/bash
#SBATCH -J s56_allele
#SBATCH -p debug
#SBATCH -n 4
#SBATCH --mem=40G
#SBATCH -t 00:59:00
#SBATCH -o /public/home/yucong/scbrain/logs/s56_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s56_%j.err
source /etc/profile
module load apps/anaconda3/2024.02
python /public/home/yucong/scbrain/scripts/s56_allele_v3.py
