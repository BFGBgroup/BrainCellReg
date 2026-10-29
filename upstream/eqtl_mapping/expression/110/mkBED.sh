#!/bin/bash
#SBATCH -N 6
#SBATCH -n 66
#SBATCH --mem=0

# 保存脑区域和细胞类型到数组
regions=("MEDU" "PUTM" "TCTX" "WHMT" "CRBL" "OCTX" "SNIG" "HIPP" "THAL" "FCTX")
cell_types=("Micro" "InN" "ExN" "Gran" "NPC" "Olig" "Endo" "OPC" "Perc" "Astro" "Purk")

# 循环遍历所有脑区域和细胞类型组合
for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
        input_file="/public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/cov/110/${region}.${cell_type}_Z.txt"
        output_file="${region}.${cell_type}.bed"

        awk 'NR==FNR {gsub(/"/,"",$1); a[$4]=$1" "$2" "$3; next} {gsub(/"/,"",$1); print a[$1], $0} END {print "#chr start end gene_id"}' /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/expression/Homo_sapiens.GRCh38.82.chr.bed "${input_file}" | sed 's/\"//g' | sed '1s/^/#chr start end gene_id /' > "${output_file}"
    done
done
