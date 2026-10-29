#!/bin/bash
#SBATCH -N 6
#SBATCH -n 66
#SBATCH --mem=0

awk '{
    if ($4 == "-") {
        printf "chr%s %d %d %s\n", $1, $3 - 1, $3, $5
    } else if ($4 == "+") {
        printf "chr%s %d %d %s\n", $1, $2 - 1, $2, $5
    }
}' /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/expression/Homo_sapiens.GRCh38.82.chr.gtf.txt > Homo_sapiens.GRCh38.82.chr.bed
