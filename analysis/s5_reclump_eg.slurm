#!/bin/bash
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --mem=110G
#SBATCH -t 59:00
#SBATCH -J scb_eg
#SBATCH -o /public/home/yucong/scbrain/logs/%x_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/%x_%j.err

# 韩教授 10-09 定：eSNP 只留 q<0.05 的 eGene 的。
# 这里：按 eGene 筛对 → 每个 SNP 取最小 P → 同参数重新去连锁（PLINK clump，250 kb，r2 0.8，1000G 欧洲）
# → 新的独立 eSNP 名单；随后用同一套脚本重做集合整理与注释抽取。non-eSNP 名单不变。
source /etc/profile
export LC_ALL=C
S=/public/home/yucong/scbrain
P=/public/home/yucong/scbrain_in/paper2
E0=$P/single-cell-eqtl/data/esnp
export SCB_OUT=$S/esnp_eg
export ESNP_DIR=$SCB_OUT/esnp
PL=$S/bin/plink
REF=$S/local_handoff/scbrain1/g1000_eur/g1000_eur
mkdir -p $ESNP_DIR/check $ESNP_DIR/clump
cd $S
date

# ---------- 0 等价性检查：用原来的 esnp2 + 我们的 PLINK 和参考，能否重现原来的 esnp3 ----------
for c in CRBL.Astro FCTX.ExN WHMT.Olig; do
  (
    $PL --bfile $REF --clump $E0/$c-esnp2.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.80 --clump-kb 250 --memory 4000 --out $ESNP_DIR/check/$c > /dev/null 2>&1
    awk 'NR>1 && $3!="" {print $3}' $ESNP_DIR/check/$c.clumped | sort > $ESNP_DIR/check/$c.mine
    grep -v '^$' $E0/$c-esnp3.txt | sort > $ESNP_DIR/check/$c.orig
    echo "CHECK $c: orig $(wc -l < $ESNP_DIR/check/$c.orig) mine $(wc -l < $ESNP_DIR/check/$c.mine) common $(comm -12 $ESNP_DIR/check/$c.orig $ESNP_DIR/check/$c.mine | wc -l)"
  ) &
done
wait
date

# ---------- 1 筛对、取最小 P、重新去连锁 ----------
regions=(CRBL FCTX HIPP MEDU OCTX PUTM SNIG TCTX THAL WHMT)
cells=(Astro Endo ExN Gran InN Micro NPC Olig OPC Perc Purk)
n=0
for r in "${regions[@]}"; do for ct in "${cells[@]}"; do
  c=$r.$ct
  (
    if [ -s $E0/$c-esnp.txt ]; then
      awk '$26 < 0.05' $E0/$c-esnp.txt > $ESNP_DIR/$c-esnp.txt
    else
      : > $ESNP_DIR/$c-esnp.txt
    fi
    awk 'BEGIN{print "SNP P"} { if (!($2 in m) || $7+0 < m[$2]+0) m[$2]=$7 } END{ for (k in m) print k, m[k] }' $ESNP_DIR/$c-esnp.txt > $ESNP_DIR/$c-esnp2.txt
    if [ $(wc -l < $ESNP_DIR/$c-esnp2.txt) -gt 1 ]; then
      $PL --bfile $REF --clump $ESNP_DIR/$c-esnp2.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.80 --clump-kb 250 --memory 4000 --out $ESNP_DIR/clump/$c > /dev/null 2>&1
      awk 'NR>1 && $3!="" {print $3}' $ESNP_DIR/clump/$c.clumped > $ESNP_DIR/$c-esnp3.txt
    else
      : > $ESNP_DIR/$c-esnp3.txt
    fi
  ) &
  n=$((n+1)); [ $((n % 18)) -eq 0 ] && wait
done; done
wait
echo "reclump done $(date)"
echo "pairs kept: $(cat $ESNP_DIR/*-esnp.txt | wc -l); esnp2 (去重前 SNP 数合计): $(cat $ESNP_DIR/*-esnp2.txt | grep -vc '^SNP P$'); esnp3 合计: $(cat $ESNP_DIR/*-esnp3.txt | wc -l); 有 esnp3 的组合: $(find $ESNP_DIR -maxdepth 1 -name '*-esnp3.txt' -size +0 | wc -l); 并集: $(cat $ESNP_DIR/*-esnp3.txt | sort -u | wc -l)"
grep -h "variants\? not found\|missing" $ESNP_DIR/clump/CRBL.Astro.log | head -3

# ---------- 2 集合整理 + 注释抽取（同一套脚本，输出到 esnp_eg） ----------
bash scripts/s2_annot.slurm
echo "S5_ALL_DONE $(date)"
