#!/bin/bash
# FAVOR chromosome 6 (missing so far: the SQL "Annotator6" file on Dataverse holds chromosome 3).
# Source: FAVOR Essential Database, CSV version, chr6.tar.gz (Dataverse file 6170516, md5 ddd75af5d441e54b774c87a9a3125afc).
# 1) keep the rows at our chr6 positions, same columns as annot/favor/chr<N>.tsv (apc_protein_function_v3 fills the
#    apc_protein_function slot; the Essential CSV has only the v3 score); 2) back up esnp_eg/results and rerun s4_stats.py;
# 3) rerun s19_ks_ctx.py for eg, lnc, pc.
#SBATCH -J scb_fav6
#SBATCH -p cyhpc_1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 03:00:00
#SBATCH -o /public/home/yucong/scbrain/logs/s23_fav6_%j.out
#SBATCH -e /public/home/yucong/scbrain/logs/s23_fav6_%j.err
set -o pipefail
source /etc/profile
export LC_ALL=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
S=/public/home/yucong/scbrain
O=$S/esnp_eg
SRC=$S/external/favor_essential_csv/chr6.tar.gz
OUT=$O/annot/favor/chr6.tsv
date
echo "md5 check"; echo "ddd75af5d441e54b774c87a9a3125afc  $SRC" | md5sum -c - || exit 1
date
tar xzf $SRC --wildcards '*.csv' -O | gawk -v POS=$O/sets/pos_chr6.txt 'BEGIN{
    FPAT = "([^,]*)|(\"[^\"]*\")"; OFS = "\t"
    while ((getline p < POS) > 0) keep[p] = 1
  }
  $1 == "variant_vcf" { next }
  ($3 in keep) {
    line = $2 OFS $3 OFS $4 OFS $5
    for (i = 6; i <= 22; i++) line = line OFS ($i == "" ? "\\N" : $i)
    line = line OFS ($26 == "" ? "\\N" : $26) OFS ($28 == "" ? "\\N" : $28) OFS ($33 == "" ? "\\N" : $33) OFS ($34 == "" ? "\\N" : $34)
    print line
  }' > $OUT.tmp || { echo EXTRACT_FAILED; exit 1; }
mv $OUT.tmp $OUT
echo "chr6 rows: $(wc -l < $OUT) | distinct positions: $(cut -f2 $OUT | sort -u | wc -l) of $(wc -l < $O/sets/pos_chr6.txt)"
date
cp -rp $O/results $O/results_before_chr6_20261015
cd $S
export SCB_OUT=$O NPROC=30
/public/home/yucong/.conda/envs/scvi/bin/python scripts/s4_stats.py || { echo S4_FAILED; exit 1; }
date
for s in eg lnc pc; do ESET=$s NPROC=16 /public/home/yucong/.conda/envs/scvi/bin/python scripts/s19_ks_ctx.py || { echo S19_FAILED; exit 1; }; done
date
echo "FAV6_DONE"
