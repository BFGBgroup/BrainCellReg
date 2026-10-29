#!/bin/bash
# Genotype boxes for every significant eSNP-gene pair (PI #118). Free debug queue (1 h cap).
#SBATCH -J scb_genoall
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 01:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s44_genoall_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s44_genoall_%j.err
PY=/public/home/yucong/.conda/envs/scvi/bin/python
S=/public/home/yucong/scbrain/scripts
mkdir -p /public/home/yucong/scbrain/tmp/duck
date
$PY $S/s44_geno_box_all.py || echo GENO_ALL_FAILED
date
