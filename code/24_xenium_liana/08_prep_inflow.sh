#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=08_prep_inflow
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o ../../processed-data/24_xenium_liana/logs/08_prep_inflow.txt
#SBATCH -e ../../processed-data/24_xenium_liana/logs/08_prep_inflow.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load liana_plus/1.7.1

## List current modules for reproducibility
module list

## Edit with your job command
python 08_prep_inflow.py

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.4.0
## available from http://research.libd.org/slurmjobs/
