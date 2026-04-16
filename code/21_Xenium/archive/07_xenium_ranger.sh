#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=128G
#SBATCH --job-name=07_xenium_ranger
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-8%8

## get sample name and data dir
read SAMPLE_NAME SAMPLE_DIR < <(awk -v i=${SLURM_ARRAY_TASK_ID} -F'\t' 'NR==i {print $1 "\t" $2}' xenium_samples.tsv)

## Explicitly pipe script output to a log
log_path=logs/07_xenium_ranger_${SAMPLE_NAME}_${SLURM_ARRAY_TASK_ID}.txt

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
module load xeniumranger/4.0.0

## List current modules for reproducibility
module list


echo "Processing sample: ${SAMPLE_NAME}"

## Edit with your job command
xeniumranger run \
  --id=${SAMPLE_NAME} \
  --xenium-bundle=/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/raw-data/xenium/${SAMPLE_DIR}\
  --output-dir=/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/21_Xenium/07_xenium_ranger/${SAMPLE_NAME} \
  --localcores=16 \
  --localmem=128


echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
