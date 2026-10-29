#!/bin/bash
#SBATCH -J scb_ldres
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 01:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s17_ldres_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s17_ldres_%j.err
# 频率判不出朝向的位点（约 94 万个）改用连锁不平衡判：s17_ld_resolve.py。
# 两道自检不过就以非零退出（后面排着的 s16a / s16 不会启动）：
#   1) 对频率已判定的位点留一法复判，与频率判定的一致率 >= 0.99
#   2) 连锁判出来的 lead 位点在 Bryois 合并单核里的方向一致率 >= 0.85
# 之后用新表重写 CRBL.Astro 的 esi，再跑一遍阿尔茨海默病 × CRBL.Astro 的 SMR，与旧 esi、只用频率的那一版对比。
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
IN=/public/home/yucong/scbrain_in
S=/public/home/yucong/scbrain/scripts
PY=/public/home/yucong/.conda/envs/scvi/bin/python
SMR=$IN/SMR_inputs/smr-1.3.1/smr
BF=/public/home/yucong/scbrain/local_handoff/scbrain1/g1000_eur/g1000_eur
GW=/public/home/yucong/scbrain/gwas10/GCST013196-smr.tsv
FIX=/public/home/yucong/scbrain/allele_fix
CX=CRBL.Astro
echo "$(date '+%F %T') resolve"
$PY $S/s17_ld_resolve.py resolve || { echo RESOLVE_FAILED; exit 1; }
echo "$(date '+%F %T') bryois"
$PY $S/s17_ld_resolve.py bryois || { echo BRYOIS_FAILED; exit 1; }
echo "$(date '+%F %T') single-context SMR"
$PY $S/s15_allele_fix.py esi $CX ld || { echo ESI_FAILED; exit 1; }
$SMR --bfile $BF --gwas-summary $GW --beqtl-summary $FIX/eqtl_ld/$CX --out $FIX/check/new_ld --diff-freq 1 --thread-num 4 > $FIX/check/new_ld.log 2>&1
grep -v -i "^warning" $FIX/check/new_ld.log | grep -i -E "included after|excluded|probes|error" | head -8
[ -s $FIX/check/new_ld.smr ] || { echo SMR_FAILED; exit 1; }
echo "== old esi vs frequency + LD"
$PY $S/s15_allele_fix.py compare $FIX/check/old.smr $FIX/check/new_ld.smr
echo "== frequency only vs frequency + LD"
$PY $S/s15_allele_fix.py compare $FIX/check/new.smr $FIX/check/new_ld.smr
echo "LDRES_ALL_DONE $(date '+%F %T')"
