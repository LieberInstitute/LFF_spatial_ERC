#!/bin/bash

#$ -cwd
#$ -l mem_free=32G,h_vmem=32G,h_fsize=100G
#$ -N BayesSpace_k_search
#$ -o ../../processed-data/05_harmony_BayesSpace/logs/03-BayesSpace_k_search.log
#$ -e ../../processed-data/05_harmony_BayesSpace/logs/03-BayesSpace_k_search.log
#$ -t 2-28
#$ -tc 15

#SBATCH -q shared
#SBATCH --mem=32G
#SBATCH --job-name=BayesSpace_k_search
#SBATCH -o ../../processed-data/05_harmony_BayesSpace/logs/03-BayesSpace_k_search.log
#SBATCH -e ../../processed-data/05_harmony_BayesSpace/logs/03-BayesSpace_k_search.log
#SBATCH --array=2-3%15

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

## List current modules for reproducibility
module load conda_R/4.3
module list

Rscript 03-BayesSpace_k_search.R

echo "**** Job ends ****"
date
