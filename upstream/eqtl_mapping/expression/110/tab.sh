#!/bin/bash
#SBATCH -N 6
#SBATCH -n 66
#SBATCH --mem=0

# 遍历所有以.bed结尾的文件
for file in /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/expression/110/*.bed; do
    # 使用sed命令将空格替换为自定义的分隔符，根据位置进行排序，将下划线改为-,并将结果保存到临时文件
    sed 's/ /\t/g' "${file}" |sort -k1,1 -k2,2n | sed 's/_/-/g' | sed '1d'|sed 's/\t\t/\t/g'> "${file}.tmp"
    # 将临时文件替换原始文件
    mv "${file}.tmp" "${file}"
done
