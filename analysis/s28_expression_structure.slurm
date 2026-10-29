#!/bin/bash
# Expression structure (s28_expression_structure.R): clustering, UMAP, CellChat on the 13,541 deconvolved profiles.
#SBATCH -J scb_exprstruct
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s28_fig2cde_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s28_fig2cde_%j.err
source /etc/profile
module load apps/R/4.5.0
cd /public/home/yucong/scbrain/scripts
date
Rscript s28_expression_structure.R || echo FIG2CDE_FAILED
date
