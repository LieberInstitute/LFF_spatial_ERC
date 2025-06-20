#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=32_sn_subcluster_hierarchical_cluster
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/32_sn_subcluster_hierarchical_cluster.txt
#SBATCH -e logs/32_sn_subcluster_hierarchical_cluster.txt
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
Rscript 32_sn_subcluster_hierarchical_cluster.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
