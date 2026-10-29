#!/bin/bash
#SBATCH -J scb_scgB
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -a 0-9
#SBATCH -o /public/home/yucong/scbrain/logs/s9b_scgB_%A_%a.out
#SBATCH -e /public/home/yucong/scbrain/logs/s9b_scgB_%A_%a.err
# scGWAS 后半段（照贺奕捷 3.<脑区>-scGWAS.sh 和 4.sh / 4.R）：一个数组任务一个脑区，
# 逐个配置文件跑 scGWAS_r1.jar（同时 K 个），再用他的 4.R 判显著模块。4.R 只删掉写死的 setwd 一行，其余不动。
IN=/public/home/yucong/scbrain_in
SG=$IN/scGWAS_inputs
HIS=$IN/paper2/eqtl-gwas/scGWAS/run
RUN=/public/home/yucong/scbrain/scgwas10
JAVA=$SG/java/jdk-25.0.1/bin/java; JAR=$SG/scGWAS-main/code/scGWAS_r1.jar
K=${K:-10}; XMX=${XMX:-8g}
REG=(MEDU PUTM TCTX WHMT CRBL OCTX SNIG HIPP THAL FCTX)
r=${REG[$SLURM_ARRAY_TASK_ID]}
source /etc/profile >/dev/null 2>&1; module load apps/R/4.5.0
[ -s $RUN/4.R ] || sed '/^setwd(/d' $HIS/4.R > $RUN/4.R.$$ && mv -n $RUN/4.R.$$ $RUN/4.R 2>/dev/null; rm -f $RUN/4.R.$$
echo "$(date '+%F %T') region $r configs $(ls $RUN/configure/$r.*-configure.txt | wc -l) K=$K XMX=$XMX"
for cfg in $(ls $RUN/configure/$r.*-configure.txt | sort); do
  b=$(basename $cfg -configure.txt); c=${b%%-GCST*}; od=$RUN/out/$c
  [ -f $od/$b.jar.done ] && continue
  ( $JAVA -Xmx$XMX -jar $JAR $cfg > $od/$b.jar.log 2>&1 && touch $od/$b.jar.done || echo "JAR_FAIL $b: $(tail -n 1 $od/$b.jar.log)" ) &
  while [ $(jobs -rp | wc -l) -ge $K ]; do sleep 2; done
done
wait
echo "$(date '+%F %T') jar done $(ls $RUN/out/$r.*/*.jar.done 2>/dev/null | wc -l)"
for z in $RUN/out/$r.*/*.combined-z.txt; do
  h=${z%.combined-z.txt}
  [ -f $h.step4.done ] && continue
  ( cd $RUN && Rscript $RUN/4.R $h $h > $h.step4.log 2>&1 && touch $h.step4.done ) &
  while [ $(jobs -rp | wc -l) -ge $K ]; do sleep 1; done
done
wait
echo "$(date '+%F %T') $r combined-z $(ls $RUN/out/$r.*/*.combined-z.txt 2>/dev/null | wc -l) sig_module files $(ls $RUN/out/$r.*/*.sig_module.txt 2>/dev/null | wc -l) nonempty $(find $RUN/out/$r.* -name '*.sig_module.txt' -size +1c | wc -l)"
echo SCGB_DONE_$r
