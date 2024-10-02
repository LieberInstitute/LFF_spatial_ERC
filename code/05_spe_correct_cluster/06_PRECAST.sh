#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=06_PRECAST
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=2-10%20

## Define loops and appropriately subset each variable for the array task ID
## all_input_genes=(HVG)
## input_genes=${all_input_genes[$(( $SLURM_ARRAY_TASK_ID / 1 % 1 ))]}
input_genes=HVG

## Explicitly pipe script output to a log
log_path=logs/06_PRECAST_${input_genes}_k${SLURM_ARRAY_TASK_ID}.txt

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
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 06_PRECAST.R --k ${SLURM_ARRAY_TASK_ID} --input_genes ${input_genes}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

