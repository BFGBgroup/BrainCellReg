#!/bin/bash
#SBATCH -J scb_esi
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 02:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s16a_esi_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s16a_esi_%j.err
# 110 个组合的 SMR .esi 按「计数等位基因」表（频率 + 连锁，allele_table_v2.tsv.gz）重写到 allele_fix/eqtl_ld，
# besd / epi 用软链接。先要有 s17 的产物。
S=/public/home/yucong/scbrain/scripts
PY=/public/home/yucong/.conda/envs/scvi/bin/python
EQ=/public/home/yucong/scbrain_in/SMR_inputs/SMR
FIX=/public/home/yucong/scbrain/allele_fix
OUT=$FIX/eqtl_ld
[ -s $FIX/allele_table_v2.tsv.gz ] || { echo NO_TABLE; exit 1; }
mkdir -p $OUT
ls $EQ/*.besd | xargs -n 1 basename | sed 's/\.besd$//' | grep -v '\.bak' > $FIX/contexts.txt
echo "contexts $(wc -l < $FIX/contexts.txt)"
one() { [ -s $OUT/$1.esi ] && [ -e $OUT/$1.esi.ld ] && [ -L $OUT/$1.besd ] && return 0; $PY $S/s15_allele_fix.py esi $1 ld || echo "ESI_FAIL $1"; }
export -f one; export PY S FIX OUT
xargs -P 6 -L 1 bash -c 'one $0' < $FIX/contexts.txt
n=$(ls $OUT/*.esi.ld 2>/dev/null | wc -l)
echo "esi written $n of $(wc -l < $FIX/contexts.txt)"
[ "$n" = "$(wc -l < $FIX/contexts.txt)" ] && echo ESI_ALL_DONE || { echo ESI_INCOMPLETE; exit 1; }
