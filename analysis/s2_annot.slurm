#!/bin/bash
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --mem=110G
#SBATCH -t 59:00
#SBATCH -J scb_annot
#SBATCH -o /public/home/yucong/scbrain/logs/%x_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/%x_%j.err

# 集合整理 + 各类注释抽取（只读 scbrain_in 和现成注释库，产物写 scbrain/esnp_func）
source /etc/profile
module load apps/bedtools/2.30.0
export LC_ALL=C
PY=/public/home/yucong/.conda/envs/scvi/bin/python
S=/public/home/yucong/scbrain
O=${SCB_OUT:-$S/esnp_func}
export SCB_OUT=$O
echo "SCB_OUT=$O ESNP_DIR=${ESNP_DIR:-default}"
P=/public/home/yucong/scbrain_in/paper2
AN=/public/home/yucong/adprog/redo/annot
FV=/public/home/yucong/scratch/favor
mkdir -p $O/sets $O/annot/enformer $O/annot/favor $O/tss $S/tmp
cd $S
date

# ---------- 1 集合 ----------
if [ ! -f $O/sets/.done ]; then
  $PY scripts/s1_sets.py && touch $O/sets/.done || { echo "s1 failed"; exit 1; }
fi
date

# 按染色体拆位置键（FAVOR 用）
awk -F'\t' -v d=$O/sets '{print $3 > (d"/pos_chr" substr($1,4) ".txt")}' $O/sets/universe_hg38.bed

# ---------- 2 TSS 距离分箱（全量 nominal） ----------
(
  [ -f $O/tss/.done ] || {
    ~/bin/duckdb -c "
      SET threads=20; SET temp_directory='$S/tmp'; SET memory_limit='40GB';
      COPY (
        SELECT regexp_extract(filename, '([A-Z]+\.[A-Za-z]+)\.cis_qtl_pairs', 1) AS ctx,
               CAST(floor(start_distance/10000.0) AS INTEGER) AS bin,
               CAST(floor(-log10(greatest(pval_nominal, 1e-300))*2) AS INTEGER) AS pbin,
               count(*) AS n
        FROM read_parquet('$P/single-cell-eqtl/data/run/new/*.parquet', filename=true)
        GROUP BY ALL
      ) TO '$O/tss/tested_pairs_bins.tsv' (DELIMITER '\t', HEADER);" && touch $O/tss/.done
  }
  echo "tss done $(date)"
) &

# ---------- 3 ENCODE cCRE ----------
(
  bedtools intersect -a $O/sets/universe_hg38.bed -b $AN/encode/GRCh38-cCREs.bed -wa -wb | cut -f4,10 > $O/annot/ccre_hits.tsv
  echo "ccre done $(date) $(wc -l < $O/annot/ccre_hits.tsv)"
) &

# ---------- 4 ReMap2022 ----------
(
  [ -f $O/annot/remap.done ] || {
    bedtools intersect -a <(zcat $AN/remap/remap2022_nr_macs2_hg38_v1_0.bed.gz) -b $O/sets/universe_hg38.bed -wa -wb \
      | awk -F'\t' -v b=$O/annot/remap_brain_hits.tsv '{
          split($4, x, ":");
          if (tolower(x[2]) ~ /brain|cortex|neuro|astro|glia|hippocamp|cerebell|microgl|oligodend/) print $13"\t"x[1]"\t"x[2] > b;
          if (!s[$13 SUBSEP x[1]]++) print $13"\t"x[1]
        }' | gzip > $O/annot/remap_hits.tsv.gz && touch $O/annot/remap.done
  }
  echo "remap done $(date)"
) &

# ---------- 5 Enformer 主成分 ----------
for c in $(seq 1 22); do
  (
    f="$AN/enformer/1000G.MAF_threshold=0.005.$c.PCA.tsv.bgz"
    zcat "$f" | awk -F'\t' 'NR==FNR{a[$1]=1; next} ($3 in a)' $O/sets/universe_rsid.txt - > $O/annot/enformer/chr$c.tsv
  ) &
done

# ---------- 6 FAVOR（文件名与染色体的对应按内容核过；6 号染色体缺） ----------
declare -A FMAP=( [1]=chr19 [2]=chr2 [3]=chr3 [4]=chr4 [5]=chr5 [7]=chr7 [8]=id6157116 [9]=chr9 [10]=id6157115 \
                  [11]=chr11 [12]=id6157112 [13]=chr12 [14]=chr13 [15]=chr14 [16]=chr15 [17]=chr16 [18]=chr17 \
                  [19]=id6157107 [20]=chr1 [21]=chr20 [22]=chr21 )
for c in "${!FMAP[@]}"; do
  (
    out=$O/annot/favor/chr$c.tsv
    [ -f $out.done ] || {
      zcat $FV/${FMAP[$c]}.sql.gz | awk -F'\t' -v C=$c 'BEGIN{OFS="\t"} NR==FNR{a[$1]=1; next}
        ($3 in a) && $2==C {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$26,$28,$33,$34}' \
        $O/sets/pos_chr$c.txt - > $out && touch $out.done
    }
    echo "favor chr$c done $(date) $(wc -l < $out)"
  ) &
done

wait
date
ls -la $O/annot $O/annot/enformer | head -60
echo ALL_DONE
