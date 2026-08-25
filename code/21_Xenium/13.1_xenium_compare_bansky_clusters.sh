#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=13.1_xenium_compare_bansky_clusters
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/13.1_xenium_compare_bansky_clusters.txt
#SBATCH -e logs/13.1_xenium_compare_bansky_clusters.txt
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
module load conda_R/4.5.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 13.1_xenium_compare_bansky_clusters.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
