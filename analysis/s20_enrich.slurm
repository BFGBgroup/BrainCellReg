#!/bin/bash
#SBATCH -J scb_enrich
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s20_enrich_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s20_enrich_%j.err
source /etc/profile
module load apps/R/4.5.0
cd /public/home/yucong/scbrain/scripts
Rscript s20_enrich.R
