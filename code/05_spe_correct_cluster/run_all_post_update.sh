#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=run_all_post_update
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/run_all_post_update.txt
#SBATCH -e logs/run_all_post_update.txt
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

## Run all (most) code post 19_SpD_update_spe

## run modeling on annotated data
rm -f "logs/20_model_pseudobulk_anno.txt"
id1=$(sbatch --parsable 20_model_pseudobulk_anno.sh)

## plot updated modeling results - dependent on 20_model_pseudobulk_anno
rm -f "logs/22_SpD_clean_plots.txt"
sbatch --parsable --dependency=afterok:$id1 22_SpD_clean_plots.sh

## pseudobulk by SpD
rm -f "logs/25_SpD_pseudobulk.txt"
id2=$(sbatch --parsable 25_SpD_pseudobulk.sh)

## mean ratio
rm -f "logs/27_SpD_MeanRatio.txt"
id3=$(sbatch --parsable 27_SpD_MeanRatio.sh)

## plot heatmaps - dependent on 20_model_pseudobulk_anno & 25_SpD_pseudobulk
rm -f "logs/23_SpD_heatmaps.txt"
sbatch --parsable --dependency=afterok:$id1:$id2:$id3 23_SpD_heatmaps.sh

echo "**** Job ends ****"
date


