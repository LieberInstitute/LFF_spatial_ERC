#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 8
#SBATCH --job-name=LC-cellranger
#SBATCH -o logs/cellranger_slurm_240528o.txt
#SBATCH --array=1-3

# Files to process
# 1c_LC_SVB 
# 3c_LC_SVB 
# 4c_LC_SVB 
 

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load CellRanger
module load cellranger/7.2.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' sample-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPLEID=$(awk 'BEGIN {FS="\t"} {print $2}' sample-list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")

echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs05/lieber/marmaypag/pilot-lc_LIBD4140/LC_pilots/raw-data/FASTQ/${SAMPLE} \
    --sample=${SAMPLEID} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs05/lieber/marmaypag/pilot-lc_LIBD4140/LC_pilots/processed-data/01_cellranger/
mv ${SAMPLE} /dcs05/lieber/marmaypag/pilot-lc_LIBD4140/LC_pilots/processed-data/01_cellranger/

echo "**** Job ends ****"
date