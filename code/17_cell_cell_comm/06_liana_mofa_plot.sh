#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=40G
#SBATCH --job-name=06_liana_mofa_plot
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o ../../processed-data/17_cell_cell_comm/logs/06_liana_mofa_plot_%a.txt
#SBATCH -e ../../processed-data/17_cell_cell_comm/logs/06_liana_mofa_plot_%a.txt
#SBATCH --array=1-3%3

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load liana_plus/1.5.1

## List current modules for reproducibility
module list

python 06_liana_mofa_plot.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
