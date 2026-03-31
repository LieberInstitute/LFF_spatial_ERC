#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=16G
#SBATCH --job-name=11_liana_DEA
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/17_cell_cell_comm/logs/11_liana_DEA.txt
#SBATCH -e ../../processed-data/17_cell_cell_comm/logs/11_liana_DEA.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

repo_dir=$(git rev-parse --show-toplevel)
in_path=$repo_dir/processed-data/13_compile_DGE/01_compile_DGE/sn_fine/DGE_results_carrier_sn_fine.Rds
out_path=$repo_dir/processed-data/13_compile_DGE/01_compile_DGE/sn_fine/DGE_results_carrier_sn_fine.csv.gz

#   Convert DEA results to CSV for input to Python
module load conda_R/4.5
Rscript -e "library(tidyverse); readRDS('$in_path') |> write_csv('$out_path')"
module unload conda_R

#   Main analysis
module load liana_plus/1.5.1
python 11_liana_DEA.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
