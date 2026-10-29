#!/bin/bash
#SBATCH -N 6
#SBATCH -n 66
#SBATCH --mem=0
# 定义处理函数
process_file() {
    full_file_name=$1
    file_name=$(echo "$full_file_name" | sed 's/_Z.txt$//')

    Rscript -e "
        gene_exp = read.table(file = '$full_file_name', header = TRUE, sep = ' ')
        gene_exp_transpose <- as.data.frame(t(gene_exp))
        write.table(gene_exp_transpose, file = paste('./output/', '$file_name', '.csv', sep = ''), sep = ',', col.names = TRUE, row.names = FALSE)
        write.table(gene_exp_transpose, file = paste('./output/', '$file_name', '.txt', sep = ''), sep = '\t', row.names = TRUE)
    "
}

export -f process_file

# 并行处理文件，每次并行处理10个文件
cat brain_cell.txt | parallel -j 10 process_file
