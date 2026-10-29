#!/bin/bash
#SBATCH -N 10
#SBATCH --ntasks 30
#SBATCH --ntasks-per-node 3
#SBATCH --cpus-per-task=2
#SBATCH --mem=0
#SBATCH --time=10-0

regions=("MEDU" "PUTM" "TCTX" "WHMT" "CRBL" "OCTX" "SNIG" "HIPP" "THAL" "FCTX")
cell_types=("Micro" "InN" "ExN" "Gran" "NPC" "Olig" "Endo" "OPC" "Perc" "Astro" "Purk")

# Function to run the Python command
run_python_command() {
    python3 -m tensorqtl /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/genotype/out /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/expression/110/${1}.${2}.bed ${1}.${2} --covariates /public/home/hanzhijie/heyijie/paper2/single-cell-eqtl/data/cov/${1}.${2}-cov.txt --mode cis
}

export -f run_python_command  # Export the function for use with GNU Parallel

# Run the loop in parallel
parallel -j 30 run_python_command ::: "${regions[@]}" ::: "${cell_types[@]}"
