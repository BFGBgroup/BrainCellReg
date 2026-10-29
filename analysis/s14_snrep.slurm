#!/bin/bash
#SBATCH -J scb_snrep
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s14_snrep_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s14_snrep_%j.err
# 单核 eQTL 复现检验（Bryois 2022）。规格：progress/分析规格_单核eQTL复现检验_v1.md
IN=/public/home/yucong/scbrain_in
S=/public/home/yucong/scbrain/scripts
PY=/public/home/yucong/.conda/envs/scvi/bin/python
PQ=$IN/paper2/single-cell-eqtl/data/run/new
BRY=/public/home/yucong/scbrain/external/bryois2022
OUT=/public/home/yucong/scbrain/sn_replication
mkdir -p $OUT/cand $OUT/bry
echo "start $(date '+%F %T'); duckdb $(~/bin/duckdb --version)"

$PY $S/s14_snrep.py prep || { echo PREP_FAILED; exit 1; }

# 候选对照：每个组合的 nominal 结果里 P>0.05 的对，按哈希固定抽约 1/300
cand() {
  pre=$1
  [ -s $OUT/cand/$pre.tsv ] && return 0
  ~/bin/duckdb -c "SET threads=2; COPY (SELECT phenotype_id AS gene, variant_id AS snp, start_distance AS dist, af, pval_nominal AS p, slope
      FROM read_parquet('$PQ/$pre.cis_qtl_pairs.chr*.parquet')
      WHERE pval_nominal > 0.05 AND (hash(phenotype_id || ':' || variant_id) % 300) = 0) TO '$OUT/cand/$pre.tsv.tmp' (DELIMITER '\t', HEADER);" > /dev/null \
    && mv $OUT/cand/$pre.tsv.tmp $OUT/cand/$pre.tsv || echo "CAND_FAIL $pre"
}
export -f cand; export OUT PQ
xargs -P 8 -L 1 bash -c 'cand $0' < $OUT/contexts.txt
echo "candidate tables $(ls $OUT/cand/*.tsv 2>/dev/null | wc -l) of $(wc -l < $OUT/contexts.txt) $(date '+%T')"

$PY $S/s14_snrep.py match || { echo MATCH_FAILED; exit 1; }

# 在 Bryois 的 9 套 × 22 条染色体里查这些对
scan() {
  f=$1; b=$(basename $f .gz)
  zcat $f | awk 'NR==FNR{k[$1" "$2]=1; next} {n=split($1,a,"_"); if ((a[n]" "$2) in k) print a[n]"\t"$2"\t"$3"\t"$4"\t"$5}' $OUT/keys.txt - > $OUT/bry/$b.tsv.tmp \
    && mv $OUT/bry/$b.tsv.tmp $OUT/bry/$b.tsv || echo "SCAN_FAIL $b"
}
export -f scan
ls $BRY/*.gz | grep -E '\.[0-9]+\.gz$' | xargs -P 16 -L 1 bash -c 'scan $0'
echo "bryois tables $(ls $OUT/bry/*.tsv 2>/dev/null | wc -l) of 198; rows $(cat $OUT/bry/*.tsv | wc -l) $(date '+%T')"

$PY $S/s14_snrep.py stats || { echo STATS_FAILED; exit 1; }
echo "SNREP_DONE $(date '+%F %T')"
