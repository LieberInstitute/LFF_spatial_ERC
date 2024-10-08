#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=12_nnSVG
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-31%20

## Define loops and appropriately subset each variable for the array task ID
all_sample=(Br1039 Br1289 Br1556 Br1691 Br1706 Br2305 Br2582 Br3974 Br5161 Br5212 Br5276 Br5367 Br5415 Br5426 Br5460 Br5517 Br5529 Br5599 Br5634 Br5712 Br5832 Br5854 Br5941 Br6085 Br6098 Br6161 Br6263 Br6321 Br6423 Br6476 Br6538)
sample=${all_sample[$(( $SLURM_ARRAY_TASK_ID / 1 % 31 ))]}

## Explicitly pipe script output to a log
log_path=logs/12_nnSVG_${sample}_${SLURM_ARRAY_TASK_ID}.txt

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
Rscript 12_nnSVG.R --sample ${sample}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

