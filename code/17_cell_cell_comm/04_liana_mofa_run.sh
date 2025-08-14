#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=64G
#SBATCH --job-name=04_liana_mofa_run
#SBATCH -c 4
#SBATCH -t 2-00:00:00
#SBATCH -o ../../processed-data/17_cell_cell_comm/logs/04_liana_mofa_run_%a.txt
#SBATCH -e ../../processed-data/17_cell_cell_comm/logs/04_liana_mofa_run_%a.txt
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

python 04_liana_mofa_run.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
