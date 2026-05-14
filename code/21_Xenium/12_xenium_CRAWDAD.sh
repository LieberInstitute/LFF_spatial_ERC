#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=12_xenium_CRAWDAD
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-16%8

## Define loops and appropriately subset each variable for the array task ID
all_brnum=(Br1039 Br1556 Br5367 Br5599 Br1706 Br5415 Br5941 Br6538 Br2582 Br5426 Br5517 Br6098 Br5161 Br5460 Br5634 Br5854)
brnum=${all_brnum[$(( $SLURM_ARRAY_TASK_ID / 1 % 16 ))]}

## Explicitly pipe script output to a log
log_path=logs/12_xenium_CRAWDAD_${brnum}_${SLURM_ARRAY_TASK_ID}.txt

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
Rscript 12_xenium_CRAWDAD.R --brnum ${brnum}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

