#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=02_optimal_cluster_voomLmFit
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-7%20

## Define loops and appropriately subset each variable for the array task ID
all_celltype=(Astro Micro OPC Oligo Vasc Excit Inhib)
celltype=${all_celltype[$(( $SLURM_ARRAY_TASK_ID / 1 % 7 ))]}

## Explicitly pipe script output to a log
log_path=logs/02_optimal_cluster_voomLmFit_${celltype}_${SLURM_ARRAY_TASK_ID}.txt

{
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
Rscript 02_optimal_cluster_voomLmFit.R --celltype ${celltype}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

