#!/bin/bash
# Rerun of part 2 only (external SMR). Job 756917 ran out of memory with 16 parallel SMR runs at 64G (171/180 runs killed).
#SBATCH -J scb_smrext
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=360G
#SBATCH -t 03:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s25_rerun_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s25_rerun_%j.err
set -o pipefail
PY=/public/home/yucong/.conda/envs/scvi/bin/python
S=/public/home/yucong/scbrain/scripts
rm -rf /public/home/yucong/scbrain/indep/smr_ext
date
/usr/bin/time -v true 2>/dev/null
NPROC=12 $PY $S/s25_smr_external.py run || echo SMR_RUN_FAILED
$PY $S/s25_smr_external.py stats || echo SMR_STATS_FAILED
date
echo SMR_EXT_DONE
