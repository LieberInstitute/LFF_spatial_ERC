#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -p katun
#SBATCH --job-name=02_spaceranger_segment
#SBATCH -o logs/02_spaceranger_segment_31.txt
#SBATCH -e logs/02_spaceranger_segment_31.txt
#SBATCH --array=17

#   Segment the full-res image for each Visium capture area

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/4.0.1

## List current modules for reproducibility
module list

repo_dir=$(git rev-parse --show-toplevel)
SAMPLE=$(cat spaceranger_parameters_1v-16v.txt spaceranger_parameters_17v-31v_trimmed.txt | awk 'BEGIN {FS=","} {print $1}' | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMG_PATH=$(cat spaceranger_parameters_1v-16v.txt spaceranger_parameters_17v-31v_trimmed.txt | awk 'BEGIN {FS=","} {print $4}' | awk "NR==${SLURM_ARRAY_TASK_ID}")
OUT_DIR=$repo_dir/processed-data/01_spaceranger/segmentation/${SAMPLE}

mkdir -p ${OUT_DIR}

echo "Processing sample ${SAMPLE}"

spaceranger segment \
    --id=${SAMPLE} \
    --tissue-image=${IMG_PATH} \
    --output-dir=${OUT_DIR} \
    --localcores=8 \
    --localmem=80 \
    --disable-ui

echo "**** Job ends ****"
date
