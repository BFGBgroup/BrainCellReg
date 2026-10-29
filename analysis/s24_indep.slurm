#!/bin/bash
# Independent replication (spec: progress/分析规格_独立数据验证_v1.md): unzip the BESD files, part 1 (GTEx region-level
# replication of the lead pairs) and part 2 (SMR with external brain eQTL).
#SBATCH -J scb_indep
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s24_indep_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s24_indep_%j.err
set -o pipefail
PY=/public/home/yucong/.conda/envs/scvi/bin/python
S=/public/home/yucong/scbrain/scripts
E=/public/home/yucong/scbrain/external/indep
cd $E
for f in *.md5sum; do z=${f%.md5sum}; sum=$(awk '{print $1}' $f); echo "$sum  $z" | md5sum -c - || echo "MD5_MISMATCH $z"; done
mkdir -p unz
for z in Brain_*.zip; do t=${z%.zip}; [ -d unz/$t ] || { mkdir -p unz/$t && unzip -q -o $z -d unz/$t; }; done
[ -d unz/BrainMeta ] || { mkdir -p unz/BrainMeta && tar xzf BrainMeta_cis_eqtl_summary.tar.gz -C unz/BrainMeta; }
find unz -name '*.besd' | sort
date
for t in Brain_Frontal_Cortex_BA9 Brain_Putamen_basal_ganglia Brain_Hippocampus Brain_Substantia_nigra Brain_Cerebellum Brain_Cortex Brain_Cerebellar_Hemisphere; do
  $PY $S/s24_indep_region.py query $t || echo "QUERY_FAILED $t"
done
$PY $S/s24_indep_region.py stats || echo STATS_FAILED
date
NPROC=16 $PY $S/s25_smr_external.py run || echo SMR_RUN_FAILED
$PY $S/s25_smr_external.py stats || echo SMR_STATS_FAILED
date
echo INDEP_DONE
