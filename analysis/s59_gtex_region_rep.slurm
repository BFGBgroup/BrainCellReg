#!/bin/bash
#SBATCH -J s59_gtex
#SBATCH -p debug
#SBATCH -n 8
#SBATCH --mem=60G
#SBATCH -t 00:59:00
#SBATCH -o /public/home/yucong/scbrain/logs/s59_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s59_%j.err
source /etc/profile
module load apps/anaconda3/2024.02
python /public/home/yucong/scbrain/scripts/s59_gtex_region_rep.py
