#!/bin/bash
#SBATCH -J scb_scgC
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 00:30:00
#SBATCH -o /public/home/yucong/scbrain/logs/s9c_scgC_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s9c_scgC_%j.err
# scGWAS 的两项核对：① 用贺奕捷的 Box-Cox 输入（阿尔茨海默病 × CRBL.Astro、帕金森病旧版 × CRBL.Astro）跑 jar，
# 对他的 combined-z.txt；② 同一份输入连跑两遍，看程序本身有没有随机性（模块、分数、4.R 判出的显著模块是否一样）。
IN=/public/home/yucong/scbrain_in
SG=$IN/scGWAS_inputs
HIS=$IN/paper2/eqtl-gwas/scGWAS/run
RUN=/public/home/yucong/scbrain/scgwas10
C=$RUN/check/repro; mkdir -p $C
JAVA=$SG/java/jdk-25.0.1/bin/java; JAR=$SG/scGWAS-main/code/scGWAS_r1.jar
source /etc/profile >/dev/null 2>&1; module load apps/R/4.5.0
for b in CRBL.Astro-GCST013196 CRBL.Astro-GCST90473307; do
  c=${b%%-GCST*}
  for rep in 1 2; do
    o=$C/$b.rep$rep
    printf 'run_model=node\nnormalization_model=calibration\nmodule_score=penalty\ngwas_node_file=%s\nscrn_expr_file=%s\nnetwork_file=%s\nr_include=0.1\nr_exclude=0.05\noutfile=%s\npermutation=true\nverbose=true\n' \
      $HIS/out/$c/$b.magma_boxcox.txt $SG/run/avg/DER-22.avg.Ensembl.tsv $SG/run/network.tsv $o > $o-configure.txt
    $JAVA -Xmx4g -jar $JAR $o-configure.txt > $o.jar.log 2>&1
    ( cd $RUN && Rscript $RUN/4.R $o $o > $o.step4.log 2>&1 )
  done
  echo "== $b"
  wc -l $HIS/out/$c/$b.combined-z.txt $C/$b.rep1.combined-z.txt $C/$b.rep2.combined-z.txt | head -3
  cmp -s $C/$b.rep1.combined-z.txt $C/$b.rep2.combined-z.txt && echo "两遍 combined-z: IDENTICAL" || echo "两遍 combined-z: DIFF"
  cmp -s $C/$b.rep1.combined-z.txt $HIS/out/$c/$b.combined-z.txt && echo "对他的 combined-z: IDENTICAL" || echo "对他的 combined-z: DIFF"
  echo "模块（第 1 列）相同的行数：rep1∩rep2 $(comm -12 <(cut -f1 $C/$b.rep1.combined-z.txt | sort) <(cut -f1 $C/$b.rep2.combined-z.txt | sort) | wc -l)；rep1∩他的 $(comm -12 <(cut -f1 $C/$b.rep1.combined-z.txt | sort) <(cut -f1 $HIS/out/$c/$b.combined-z.txt | sort) | wc -l)"
  echo "前 6 列（模块、种子、边、分数）相同的行数：rep1∩他的 $(comm -12 <(cut -f1-6 $C/$b.rep1.combined-z.txt | sort) <(cut -f1-6 $HIS/out/$c/$b.combined-z.txt | sort) | wc -l)"
  for x in rep1 rep2; do echo "sig_module $x: $(( $(wc -l < $C/$b.$x.sig_module.txt 2>/dev/null || echo 1) - 1 )) 行"; done
  echo "sig_module 他的: $(( $(wc -l < $HIS/out/$c/$b.sig_module.txt) - 1 )) 行"
  echo "显著模块相同：rep1∩rep2 $(comm -12 <(cut -f1 $C/$b.rep1.sig_module.txt | sort) <(cut -f1 $C/$b.rep2.sig_module.txt | sort) | grep -vc module_genes)；rep1∩他的 $(comm -12 <(cut -f1 $C/$b.rep1.sig_module.txt | sort) <(cut -f1 $HIS/out/$c/$b.sig_module.txt | sort) | grep -vc module_genes)"
done
echo SCGC_DONE
