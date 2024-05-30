#!/bin/bash
#SBATCH --mem=60G
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg/code/logs/countNuclei_untrimmed_%a.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/code/VistoSeg/code/logs/countNuclei_untrimmed_%a.txt
#SBATCH --array=11
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
#samplelist="countNuclei_list_untrimmed.txt"

## Read inputs from countNuclei_list.txt file
FILE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" countNuclei_list_untrimmed.txt)
mask=$(echo ${FILE} | cut -d "," -f 1)
echo "using nuclei mat file ${mask}"
jsonname=$(echo ${FILE} | cut -d "," -f 2)
echo "using json file ${jsonname}"
posname=$(echo ${FILE} | cut -d "," -f 3)
echo " and using tissueposition file ${posname}"


## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname','$posname')"

echo "**** Job ends ****"
date