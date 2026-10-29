#!/bin/bash
#SBATCH -J s57_tab3
#SBATCH -p debug
#SBATCH -n 4
#SBATCH --mem=60G
#SBATCH -t 00:59:00
#SBATCH -o /public/home/yucong/scbrain/logs/s57_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s57_%j.err
source /etc/profile
module load apps/anaconda3/2024.02
python /public/home/yucong/scbrain/scripts/s57_tables_allele_v3.py
