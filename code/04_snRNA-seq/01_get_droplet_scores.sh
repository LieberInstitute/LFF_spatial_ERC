#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=01_get_droplet_scores
#SBATCH -c 1
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-31%20

## Define loops and appropriately subset each variable for the array task ID
all_sample_id=(10c_ERC_SVB 11c_ERC_SVB 12c_ERC_SVB 13c_ERC_SVB 14c_ERC_SVB 15c_ERC_SVB 16c_ERC_SVB 17c_ERC_SVB 18c_ERC_SVB 19c_ERC_SVB 1c_ERC_SVB 20c_ERC_SVB 21c_ERC_SVB 22c_ERC_SVB 23c_ERC_SVB 24c_ERC_SVB 25c_ERC_SVB 26c_ERC_SVB 27c_ERC_SVB 28c_ERC_SVB 29c_ERC_SVB 2c_ERC_SVB 30c_ERC_SVB 31c_ERC_SVB 3c_ERC_SVB 4c_ERC_SVB 5c_ERC_SVB 6c_ERC_SVB 7c_ERC_SVB 8c_ERC_SVB 9c_ERC_SVB)
sample_id=${all_sample_id[$(( $SLURM_ARRAY_TASK_ID / 1 % 31 ))]}

## Explicitly pipe script output to a log
log_path=logs/01_get_droplet_scores_${sample_id}_${SLURM_ARRAY_TASK_ID}.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_get_droplet_scores.R --sample_id ${sample_id}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/

