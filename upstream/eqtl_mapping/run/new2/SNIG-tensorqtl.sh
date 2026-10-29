#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks 11
#SBATCH --mem=0
#SBATCH --time=10-0

module load apps/R/4.3.1
cell_types=("Micro" "InN" "ExN" "Gran" "NPC" "Olig" "Endo" "OPC" "Perc" "Astro" "Purk")

# Function to run the Python command
run_python_command() {
    python3 -m tensorqtl /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/genotype/out /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/expression/110/SNIG.${1}.bed SNIG.${1} --covariates /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/cov/SNIG.${1}-cov.txt --mode cis
}

export -f run_python_command  # Export the function for use with GNU Parallel

# Run the loop in parallel
parallel -j 11 run_python_command ::: "${regions[@]}" ::: "${cell_types[@]}"
