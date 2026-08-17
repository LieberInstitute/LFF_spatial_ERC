#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=05.5_qTune_seedTest
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --array=1-10
#SBATCH -o logs/05.5_qTune_seedTest_%a.txt
#SBATCH -e logs/05.5_qTune_seedTest_%a.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Array task id: ${SLURM_ARRAY_TASK_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

## Load the R module
module load conda_R/4.5.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 05.5_qTune_seedTest.R --spe spe_PCA_SVG --dimred HARMONY --name SVGm

echo "**** Job ends ****"
date

