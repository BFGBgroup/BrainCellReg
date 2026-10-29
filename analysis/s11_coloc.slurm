#!/bin/bash
#SBATCH -J scb_coloc
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=100G
#SBATCH -t 06:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s11_coloc_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s11_coloc_%j.err
# 共定位。MODE=check：只用没换 GWAS 的 3 种病、拿贺奕捷的 smr-sig-gene.txt 核对能否重现他的结果；
# MODE=full：先汇总 SMR（s10），再对全部 q_SMR<0.05 的（组合、GWAS、基因）做共定位。
MODE=${MODE:-full}
IN=/public/home/yucong/scbrain_in
S=/public/home/yucong/scbrain/scripts
PY=/public/home/yucong/.conda/envs/scvi/bin/python
PQ=$IN/paper2/single-cell-eqtl/data/run/new
source /etc/profile >/dev/null 2>&1; module load apps/R/4.5.0
if [ "$MODE" = check ]; then
  OUT=/public/home/yucong/scbrain/coloc10_check; mkdir -p $OUT/eqtl
  grep -E 'GCST013196|GCST90014122|GCST90477562' $IN/paper2/eqtl-gwas/SMR/10disease-SMR/smr-sig-gene.txt > $OUT/smr-sig-gene.txt
else
  OUT=/public/home/yucong/scbrain/coloc10; mkdir -p $OUT/eqtl
  $PY $S/s10_smr_post.py || { echo SMR_POST_FAILED; exit 1; }
  cp /public/home/yucong/scbrain/smr10/smr-sig-gene.txt $OUT/smr-sig-gene.txt
fi
SIG=$OUT/smr-sig-gene.txt
echo "triplets $(wc -l < $SIG); contexts $(cut -f1 $SIG | sort -u | wc -l); context-gene pairs $(cut -f1,3 $SIG | sort -u | wc -l)"

# 每个组合一次性从 nominal parquet 里取出要用的基因的全部 SNP
one_ctx() {
  cx=$1; pre=${cx/-/.}
  [ -s $OUT/eqtl/$cx.tsv ] && return 0
  genes=$(awk -F'\t' -v c=$cx '$1==c{print $3}' $SIG | sort -u | sed "s/.*/'&'/" | paste -sd, -)
  ~/bin/duckdb -c "SET threads=2; COPY (SELECT phenotype_id AS genes, variant_id AS snps, start_distance AS tss, pval_nominal AS pvalue, slope AS beta
      FROM read_parquet('$PQ/$pre.cis_qtl_pairs.chr*.parquet') WHERE phenotype_id IN ($genes) AND pval_nominal IS NOT NULL AND NOT isnan(pval_nominal)) TO '$OUT/eqtl/$cx.tsv.tmp' (DELIMITER '\t', HEADER);" > /dev/null \
    && mv $OUT/eqtl/$cx.tsv.tmp $OUT/eqtl/$cx.tsv || echo "EXTRACT_FAIL $cx"
}
export -f one_ctx; export OUT SIG PQ
cut -f1 $SIG | sort -u | xargs -P 6 -L 1 bash -c 'one_ctx $0'
echo "eqtl tables $(ls $OUT/eqtl/*.tsv 2>/dev/null | wc -l) $(date '+%T')"

Rscript $S/s11_coloc.R $SIG $OUT/eqtl $OUT 10
echo "COLOC_JOB_DONE $MODE"
