#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg/code/logs/refineVNS_%a.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg/code/logs/refineVNS_%a.txt
#SBATCH --array=25-31%8
#SBATCH --mail-user=heenadivecha@gmail.com
 
echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab/R2023a

## Load toolbox for VistoSeg
toolbox='/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg'
#samplelist="refineVNS.txt"

## Read inputs from refineVNS.txt file
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" refineVNS.txt)
fname=$(echo ${SAMPLE} | cut -d "," -f 1)
echo "Processing image ${fname}"
M=$(echo ${SAMPLE} | cut -d "," -f 2)
echo "Cluster number ${M}"

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), refineVNS('$fname',$M)"

echo "**** Job ends ****"
date