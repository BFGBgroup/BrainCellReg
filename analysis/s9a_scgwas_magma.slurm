#!/bin/bash
#SBATCH -J scb_scgA
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s9a_scgA_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s9a_scgA_%j.err
# scGWAS 前半段（照贺奕捷 genes_annot/1.sh、2.sh，run/1.sh、2.sh、111.generate_all_configure.sh）：
#   基因—eSNP 注释表（只留 eGene 的 eSNP）→ MAGMA → Box-Cox → 配置文件；最后单跑一个 scGWAS 量时间。
# 先做两项核对：① 用他的旧注释表 + 阿尔茨海默病 GWAS 重现他的 MAGMA 结果；② 汇总数据只留注释表里的 SNP 与用全表结果一致。
IN=/public/home/yucong/scbrain_in
SG=$IN/scGWAS_inputs
HIS=$IN/paper2/eqtl-gwas/scGWAS/run
MAGMA=$SG/magma/magma
BF=/public/home/yucong/scbrain/local_handoff/scbrain1/g1000_eur/g1000_eur
GW=/public/home/yucong/scbrain/gwas10
ES=/public/home/yucong/scbrain/esnp_eg/esnp
GTF=$IN/paper2/single-cell-eqtl/data/expression/Homo_sapiens.GRCh38.82.chr.gtf.txt
RUN=/public/home/yucong/scbrain/scgwas10
NP=${SLURM_CPUS_PER_TASK:-24}
source /etc/profile >/dev/null 2>&1; module load apps/R/4.5.0
mkdir -p $RUN/genes_annot $RUN/GWAS $RUN/out $RUN/check $RUN/configure
export MAGMA BF RUN GW HIS

echo "== [1] 注释表 $(date '+%T')"
for f in $ES/*-esnp.txt; do
  ctx=$(basename $f -esnp.txt)
  [ -s $f ] || { echo "empty $ctx"; continue; }
  awk '{arr[$1] = arr[$1] "\t" $2} END{for (g in arr) print g arr[g]}' $f > $RUN/genes_annot/$ctx-esnp.genes.annot.txt
  awk -F'\t' 'NR==FNR {coord[$5] = $1 ":" $2 ":" $3; next}
       {s = ""; for (i = 2; i <= NF; i++) if ($i != "") s = (s == "" ? $i : s OFS $i); print $1, (($1 in coord) ? coord[$1] : "NA"), s}' OFS='\t' \
       $GTF $RUN/genes_annot/$ctx-esnp.genes.annot.txt > $RUN/genes_annot/$ctx-esnp.genes.annot.with.pos.txt
done
echo "注释表个数: $(ls $RUN/genes_annot/*.with.pos.txt | wc -l)；基因行合计: $(cat $RUN/genes_annot/*.with.pos.txt | wc -l)；没有坐标的基因行: $(cat $RUN/genes_annot/*.with.pos.txt | awk -F'\t' '$2=="NA"' | wc -l)"
echo "他的旧注释表基因行合计: $(cat $HIS/genes_annot/*.with.pos.txt | wc -l)"

echo "== [2] 汇总数据只留注释表里的 SNP $(date '+%T')"
cat $RUN/genes_annot/*.with.pos.txt $HIS/genes_annot/CRBL.Astro-esnp.genes.annot.with.pos.txt | cut -f3- | tr '\t' '\n' | sort -u > $RUN/GWAS/annot_snps.txt
wc -l < $RUN/GWAS/annot_snps.txt
cp $GW/gwasinfo.txt $RUN/GWAS/gwasinfo.txt
while read -r g N; do
  ( awk 'NR==FNR {k[$1]; next} FNR==1 || ($1 in k)' $RUN/GWAS/annot_snps.txt $GW/$g.sumstats.txt > $RUN/GWAS/$g.sumstats.txt; echo "$g N=$N kept $(($(wc -l < $RUN/GWAS/$g.sumstats.txt) - 1))" ) &
done < $RUN/GWAS/gwasinfo.txt
wait

echo "== [3] 两项核对 $(date '+%T')"
c=CRBL.Astro; g=GCST013196; N=$(awk -v g=$g '$1==g{print $2}' $RUN/GWAS/gwasinfo.txt)
$MAGMA --bfile $BF --pval $GW/$g.sumstats.txt N=$N --gene-annot $HIS/genes_annot/$c-esnp.genes.annot.with.pos.txt --out $RUN/check/his_full > $RUN/check/his_full.log 2>&1 &
$MAGMA --bfile $BF --pval $RUN/GWAS/$g.sumstats.txt N=$N --gene-annot $HIS/genes_annot/$c-esnp.genes.annot.with.pos.txt --out $RUN/check/his_filt > $RUN/check/his_filt.log 2>&1 &
wait
cmp -s $RUN/check/his_full.genes.out $HIS/out/$c/$c-$g.Magma.genes.out && echo "CHECK1 重现他的 MAGMA 结果: IDENTICAL" || { echo "CHECK1 DIFF"; diff $RUN/check/his_full.genes.out $HIS/out/$c/$c-$g.Magma.genes.out | head -6; }
cmp -s $RUN/check/his_full.genes.out $RUN/check/his_filt.genes.out && echo "CHECK2 全表与只留注释 SNP: IDENTICAL" || { echo "CHECK2 DIFF"; diff $RUN/check/his_full.genes.out $RUN/check/his_filt.genes.out | head -6; }
if ! cmp -s $RUN/check/his_full.genes.out $RUN/check/his_filt.genes.out; then
  echo "只留注释 SNP 的结果与全表不同，改回用全表"; for g2 in $(cut -f1 $RUN/GWAS/gwasinfo.txt); do ln -sf $GW/$g2.sumstats.txt $RUN/GWAS/$g2.sumstats.txt; done
fi

echo "== [4] MAGMA 全部 $(date '+%T')"
one_magma() {
  c=$1; g=$2; N=$3; od=$RUN/out/$c; mkdir -p $od
  [ -s $od/$c-$g.Magma.genes.out ] && return 0
  $MAGMA --bfile $BF --pval $RUN/GWAS/$g.sumstats.txt N=$N --gene-annot $RUN/genes_annot/$c-esnp.genes.annot.with.pos.txt --out $od/$c-$g.Magma > $od/$c-$g.Magma.runlog 2>&1
}
export -f one_magma
for a in $RUN/genes_annot/*-esnp.genes.annot.with.pos.txt; do
  c=$(basename $a -esnp.genes.annot.with.pos.txt)
  while read -r g N; do echo "$c $g $N"; done < $RUN/GWAS/gwasinfo.txt
done > $RUN/magma_tasks.txt
wc -l < $RUN/magma_tasks.txt
xargs -P $NP -L 1 bash -c 'one_magma $0 $1 $2' < $RUN/magma_tasks.txt
echo "genes.out 个数: $(ls $RUN/out/*/*.Magma.genes.out 2>/dev/null | wc -l)"

echo "== [5] Box-Cox $(date '+%T')"
one_box() { f=$1; o=${f%.Magma.genes.out}; [ -s $o.magma_boxcox.txt ] || Rscript $HIS/2.magma_boxcox.R $f $o.magma_boxcox.txt $o.magma_boxcox.pdf > $o.boxcox.runlog 2>&1; }
export -f one_box
ls $RUN/out/*/*.Magma.genes.out | xargs -P $NP -L 1 bash -c 'one_box $0'
echo "boxcox 个数: $(ls $RUN/out/*/*.magma_boxcox.txt 2>/dev/null | wc -l)"

echo "== [6] 配置文件 $(date '+%T')"
for b in $RUN/out/*/*.magma_boxcox.txt; do
  o=${b%.magma_boxcox.txt}
  printf 'run_model=node\nnormalization_model=calibration\nmodule_score=penalty\ngwas_node_file=%s\nscrn_expr_file=%s\nnetwork_file=%s\nr_include=0.1\nr_exclude=0.05\noutfile=%s\npermutation=true\nverbose=true\n' \
    $b $SG/run/avg/DER-22.avg.Ensembl.tsv $SG/run/network.tsv $o > $RUN/configure/$(basename $o)-configure.txt
done
echo "配置文件个数: $(ls $RUN/configure/*-configure.txt | wc -l)"

echo "== [7] 单跑一个 scGWAS 量时间 $(date '+%T')"
JAVA=$SG/java/jdk-25.0.1/bin/java; JAR=$SG/scGWAS-main/code/scGWAS_r1.jar
/usr/bin/time -v $JAVA -jar $JAR $RUN/configure/CRBL.Astro-GCST013196-configure.txt > $RUN/check/java_test.log 2> $RUN/check/java_test.time
echo "java exit $?"; tail -n 5 $RUN/check/java_test.log
grep -E "Elapsed|Maximum resident|Percent of CPU" $RUN/check/java_test.time
ls -la $RUN/out/CRBL.Astro/ | grep GCST013196 | cut -c25-
echo SCGA_DONE
