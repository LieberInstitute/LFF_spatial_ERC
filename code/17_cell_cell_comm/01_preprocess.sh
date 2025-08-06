#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=01_preprocess
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/17_cell_cell_comm/logs/01_preprocess_%a.txt
#SBATCH -e ../../processed-data/17_cell_cell_comm/logs/01_preprocess_%a.txt
#SBATCH --array=1-31%10

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0

## List current modules for reproducibility
module list

repo_dir=$(git rev-parse --show-toplevel)
sample_id_path_1=$repo_dir/code/01_spaceranger/spaceranger_parameters_1v-16v.txt
sample_id_path_2=$repo_dir/code/01_spaceranger/spaceranger_parameters_17v-31v_trimmed.txt
database_path=/jhpce/shared/libd/core/visium_hd/1.0/CellNEST/database/CellNEST_database.csv

sample_id=$(cat ${sample_id_path_1} ${sample_id_path_2} | awk 'BEGIN {FS=","} {print $1}' | awk "NR==${SLURM_ARRAY_TASK_ID}")
in_dir=$repo_dir/processed-data/01_spaceranger/$sample_id/outs
out_dir=$repo_dir/processed-data/17_cell_cell_comm

cellnest preprocess \
    --data_name=${sample_id} \
    --data_from=${in_dir} \
    --data_to=${out_dir}/input_graph/${sample_id}/ \
    --metadata_to=${out_dir}/metadata/${sample_id}/ \
    --database_path=${database_path}

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
