#!/bin/bash
#SBATCH -J scb_fixdn
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=100G
#SBATCH -t 03:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s18_fixdn_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s18_fixdn_%j.err
# 改正等位基因后的 SMR（smr10_fix）往下：汇总并与原结果对比 → 给新进 q_SMR<0.05 的三元组补共定位
# （共定位不用等位基因字母，已有的结果直接沿用）→ 三路并表到 integrate10_fix → 与 integrate10 对比。
# 原来的 smr10 / coloc10 / integrate10 都不动。
IN=/public/home/yucong/scbrain_in
B=/public/home/yucong/scbrain
S=$B/scripts
PY=/public/home/yucong/.conda/envs/scvi/bin/python
PQ=$IN/paper2/single-cell-eqtl/data/run/new
source /etc/profile >/dev/null 2>&1; module load apps/R/4.5.0
n=$(ls $B/smr10_fix/*/*.done 2>/dev/null | wc -l)
echo "SMR results with .done: $n (expected 1100 minus contexts without data)"
$PY $S/s18_smr_fix_post.py || { echo SMR_FIX_POST_FAILED; exit 1; }

OUT=$B/coloc10_fix; mkdir -p $OUT/eqtl
SIG=$OUT/smr-sig-gene.txt
if [ -s $SIG ]; then
  echo "coloc to run: triplets $(wc -l < $SIG); contexts $(cut -f1 $SIG | sort -u | wc -l)"
  one_ctx() {
    cx=$1; pre=${cx/-/.}
    genes=$(awk -F'\t' -v c=$cx '$1==c{print $3}' $SIG | sort -u | sed "s/.*/'&'/" | paste -sd, -)
    ~/bin/duckdb -c "SET threads=2; COPY (SELECT phenotype_id AS genes, variant_id AS snps, start_distance AS tss, pval_nominal AS pvalue, slope AS beta
        FROM read_parquet('$PQ/$pre.cis_qtl_pairs.chr*.parquet') WHERE phenotype_id IN ($genes) AND pval_nominal IS NOT NULL AND NOT isnan(pval_nominal)) TO '$OUT/eqtl/$cx.tsv.tmp' (DELIMITER '\t', HEADER);" > /dev/null \
      && mv $OUT/eqtl/$cx.tsv.tmp $OUT/eqtl/$cx.tsv || echo "EXTRACT_FAIL $cx"
  }
  export -f one_ctx; export OUT SIG PQ
  cut -f1 $SIG | sort -u | xargs -P 6 -L 1 bash -c 'one_ctx $0'
  Rscript $S/s11_coloc.R $SIG $OUT/eqtl $OUT 10 || { echo COLOC_FAILED; exit 1; }
  mv $OUT/coloc_all.tsv $OUT/coloc_new_only.tsv
else
  echo "no new triplets need coloc"; rm -f $OUT/coloc_new_only.tsv
fi
$PY $S/s18b_fix_merge.py merge || { echo MERGE_FAILED; exit 1; }
SMR_DIR=$B/smr10_fix COLOC_DIR=$B/coloc10_fix INT_OUT=$B/integrate10_fix $PY $S/s13_integrate.py || { echo INTEGRATE_FAILED; exit 1; }
$PY $S/s18b_fix_merge.py compare
echo "FIX_DOWNSTREAM_DONE $(date '+%F %T')"
