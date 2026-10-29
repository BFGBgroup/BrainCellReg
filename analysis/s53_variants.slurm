#!/bin/bash
#SBATCH -J s53_var
#SBATCH -p debug
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH -o /public/home/yucong/scbrain/logs/s53_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s53_%j.err
source /etc/profile
module load apps/anaconda3/2024.02
python /public/home/yucong/scbrain/scripts/s53_variants.py
