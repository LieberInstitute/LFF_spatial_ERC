#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=06_Clusterwise_voomLmFit_Xenium_cell_type_anno_SpX
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/06_Clusterwise_voomLmFit_Xenium_cell_type_anno_SpX.txt
#SBATCH -e logs/06_Clusterwise_voomLmFit_Xenium_cell_type_anno_SpX.txt
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
Rscript 06_Clusterwise_voomLmFit_Xenium.R --datatype Xenium_cell_type_anno_SpX

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
