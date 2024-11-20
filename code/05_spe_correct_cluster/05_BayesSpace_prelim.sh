#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=05_BayesSpace_prelim
#SBATCH -c 1
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=2-10%10

## Explicitly pipe script output to a log
log_path=logs/05_BayesSpace_prelim_k${SLURM_ARRAY_TASK_ID}.txt

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
module load conda_R/4.4

## List current modules for reproducibility
module list

## Run BayesSpace on regular (Harmony corrected) PCAs
Rscript 05_BayesSpace_prelim.R --spe spe_postQC --dimred HARMONY --name PCA_Harmony_prelim --prelim 10x
# Rscript 05_BayesSpace_prelim.R --spe spe_postQC --dimred HARMONY --name PCA_Harmony_prelim_id --prelim 10x_id

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/
