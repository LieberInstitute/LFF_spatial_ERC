#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg/code/logs/VNS_%a.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg/code/logs/VNS_%a.txt
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
samplelist="VNS.txt"

## Read inputs from VNS.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")


## Run VNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"

echo "**** Job ends ****"
date