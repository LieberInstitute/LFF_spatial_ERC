#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=200G
#SBATCH --job-name=09_xenium_label_transfer_RCDT
#SBATCH -c 1
#SBATCH -t 3-00:00:00
#SBATCH -o logs/09_xenium_label_transfer_RCDT.txt
#SBATCH -e logs/09_xenium_label_transfer_RCDT.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 09_xenium_label_transfer_RCDT.R --cell_type_col cell_type_anno

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
