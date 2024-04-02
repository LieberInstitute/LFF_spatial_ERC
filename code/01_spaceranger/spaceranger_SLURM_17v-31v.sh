#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 8
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/01_spaceranger/logs/spaceranger_%a.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/01_spaceranger/logs/spaceranger_%a.txt
#SBATCH --array=1
#SBATCH --mail-user=heenadivecha@gmail.com
 
echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/2.1.0

## List current modules for reproducibility
module list

## Read inputs from spaceranger_parameters.txt file
FILE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" spaceranger_parameters.txt)
SAMPLE=$(echo ${FILE} | cut -d "," -f 1)
SLIDE=$(echo ${FILE} | cut -d "," -f 2)
CAPTUREAREA=$(echo ${FILE} | cut -d "," -f 3)
IMAGEPATH=$(echo ${FILE} | cut -d "," -f 4)
LOUPEPATH=$(echo ${FILE} | cut -d "," -f 5)
FASTQPATH=$(echo ${FILE} | cut -d "," -f 6)

echo "Processing sample ${SAMPLE} from slide ${SLIDE} and capture area ${CAPTUREAREA} with image ${IMAGEPATH} and aligned with ${LOUPEPATH} with FASTQs: ${FASTQPATH}"
date

## For keeping track of dates of the input files
ls -lh ${IMAGEPATH}
ls -lh ${LOUPEPATH}



## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

## Removed '--unknown-slide=visium-1' argument as the loupe files were created without slide information but Leo helped resave the file with slide information. Leo created 01_update_json_files.R for the same. 
## Added '--r1-length=26' to fix length variations in R1 FASTQ file. 
## Added '--loupe-alignment=${LOUPEPATH}' to as we now have the slide information

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --image=${IMAGEPATH} \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=${LOUPEPATH} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64 \
    --r1-length=26

## Move output
echo "Moving results to new location"
date
mkdir -p /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/01_spaceranger/
mv ${SAMPLE} /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/01_spaceranger/

echo "**** Job ends ****"
date

