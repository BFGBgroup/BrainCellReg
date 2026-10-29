#!/bin/bash
#SBATCH -J scb_smr
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=170G
#SBATCH -t 08:00:00
#SBATCH -a 0-6
#SBATCH -o /public/home/yucong/scbrain/logs/s8_smr_%A_%a.out
#SBATCH -e /public/home/yucong/scbrain/logs/s8_smr_%A_%a.err
# 新换的 7 份 GWAS × 110 个组合的 SMR，参数照贺奕捷 2-<脑区>-smr.sh：--diff-freq 1，其余默认。
# 一个数组任务一份 GWAS；同时跑 K 个组合；已有 .done 标记的跳过（可断点续跑）。
IN=/public/home/yucong/scbrain_in
SMR=$IN/SMR_inputs/smr-1.3.1/smr
BF=/public/home/yucong/scbrain/local_handoff/scbrain1/g1000_eur/g1000_eur
EQ=$IN/SMR_inputs/SMR
GW=/public/home/yucong/scbrain/gwas10
OUT=/public/home/yucong/scbrain/smr10
K=${K:-8}; T=${T:-1}
IDS=(GCST90432156 GCST003566 GCST90726344 GCST009325 GCST90018844 GCST90128471 GCST011179)
g=${IDS[$SLURM_ARRAY_TASK_ID]}
[ -s $GW/$g-smr.tsv ] || { echo "NO_GWAS $g"; exit 1; }
echo "$(date '+%F %T') start $g K=$K T=$T"
for r in MEDU PUTM TCTX WHMT CRBL OCTX SNIG HIPP THAL FCTX; do
  for c in Endo Purk Perc NPC OPC InN ExN Gran Olig Astro Micro; do
    [ -f $EQ/$r.$c.besd ] || { echo "skip missing $r.$c"; continue; }
    od=$OUT/$r-$c; mkdir -p $od
    [ -f $od/$g.done ] && continue
    (
      $SMR --bfile $BF --gwas-summary $GW/$g-smr.tsv --beqtl-summary $EQ/$r.$c --out $od/$g --diff-freq 1 --thread-num $T > $od/$g.log 2>&1
      rc=$?
      if [ $rc -eq 0 ] && [ -s $od/$g.smr ]; then touch $od/$g.done; else echo "FAIL $r.$c rc=$rc: $(tail -n 2 $od/$g.log | tr '\n' ' ')"; fi
    ) &
    while [ $(jobs -rp | wc -l) -ge $K ]; do sleep 2; done
  done
done
wait
n=$(ls $OUT/*/$g.done 2>/dev/null | wc -l)
echo "$(date '+%F %T') $g done $n / 110"
[ "$n" = 110 ] && echo SMR_DONE_$g || echo SMR_INCOMPLETE_$g
