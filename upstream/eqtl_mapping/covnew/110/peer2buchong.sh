#!/bin/bash
#SBATCH -N 6
#SBATCH -n 66
#SBATCH --mem=0
# 定义处理函数
process_file() {
    full_file_name=$1
    file_name=$(echo "$full_file_name" | sed 's/_Z.txt$//')
    peertool -f ./output/$file_name.csv -n 15 --has_header -o ./output/$file_name
}
export -f process_file
# 并行处理文件，每次并行处理10个文件
cat brain_cellbuchong.txt | parallel -j 10 process_file
