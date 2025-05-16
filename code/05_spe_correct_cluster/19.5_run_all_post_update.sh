#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=19_SpD_update_spe
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/19_SpD_update_spe.txt
#SBATCH -e logs/19_SpD_update_spe.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

MAINDIR=$(git rev-parse --show-toplevel)
CODEDIR="${MAINDIR}/code/05_spe_correct_cluster/"

cd ${CODEDIR}

## run modeling on annotated data
rm -f "logs/20_model_pseudobulk_anno.txt"
id1=$(sbatch --parsable 20_model_pseudobulk_anno.sh)

## plot updated modeling results - dependent on 20_model_pseudobulk_anno
rm -f "logs/22_SpD_clean_plots.txt"
sbatch --parsable --dependency=afterok:$id1 22_SpD_clean_plots.sh

## pseudobulk by SpD
id2=$(sbatch --parsable 25_SpD_pseudobulk.sh)

## mean ratio
sbatch --parsable 27_SpD_MeanRatio.sh

## plot heatmaps - dependent on 20_model_pseudobulk_anno & 25_SpD_pseudobulk
sbatch --parsable --dependency=afterok:$id1 23_SpD_heatmaps.sh



28_SpD_pseudobulk_pca.R

