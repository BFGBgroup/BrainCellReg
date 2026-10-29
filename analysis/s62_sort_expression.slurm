#!/bin/bash
#SBATCH -J s62_sortexpr
#SBATCH -p debug
#SBATCH -n 8
#SBATCH --mem=60G
#SBATCH -t 00:59:00
#SBATCH -o /public/home/yucong/scbrain/logs/s62_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s62_%j.err
source /etc/profile
module load apps/anaconda3/2024.02
python /public/home/yucong/scbrain/scripts/s62_sort_expression.py
