#!/bin/bash
# Database step 2: tables for the draw-on-demand plots (s26_db_plot_tables.py).
#SBATCH -J scb_dbplot
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s26_dbplot_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s26_dbplot_%j.err
PY=/public/home/yucong/.conda/envs/scvi/bin/python
S=/public/home/yucong/scbrain/scripts
mkdir -p /public/home/yucong/scbrain/tmp/duck
rm -f /public/home/yucong/scbrain/db/tables/CHECK2.txt
date
$PY $S/s26_db_plot_tables.py demo GPNMB,ADAM15,CD38,CLU || echo DEMO_FAILED
date
$PY $S/s26_db_plot_tables.py geno || echo GENO_FAILED
date
$PY $S/s26_db_plot_tables.py expr || echo EXPR_FAILED
date
echo S26_DONE
