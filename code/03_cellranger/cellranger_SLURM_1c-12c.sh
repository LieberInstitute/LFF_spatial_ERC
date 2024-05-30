#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 8
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/01_spaceranger/logs/cellranger_%a.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/01_spaceranger/logs/cellranger_%a.txt
#SBATCH --array=1-12
#SBATCH --mail-user=heenadivecha@gmail.com

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

echo "**** Job starts ****"
date

## load CellRanger
module load cellranger/7.2.0

## List current modules for reproducibility
module list

## Read inputs from cellranger_SLURM_1c-12c.txt file
FILE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" cellranger_SLURM_1c-12c.txt)
SAMPLE=$(echo ${FILE} | cut -d "," -f 1)
SAMPLEID=$(echo ${FILE} | cut -d "," -f 2)
FASTQPATH=$(echo ${FILE} | cut -d "," -f 3)


echo "Processing sample ${SAMPLE} with sampleid: ${SAMPLEID} and FASTQs: ${FASTQPATH}"
date

## For keeping track of dates of the input files
ls -lh ${SAMPLE}
ls -lh ${FASTQPATH}

## Run CellRanger
cellranger count 
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --sampleid=${SAMPLEID} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64 \

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/03_cellranger/
mv ${SAMPLE} /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/03_cellranger/

echo "**** Job ends ****"
date