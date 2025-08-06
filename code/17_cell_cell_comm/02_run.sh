#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=02_run_cpu
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/17_cell_cell_comm/logs/02_run_cpu_%a.txt
#SBATCH -e ../../processed-data/17_cell_cell_comm/logs/02_run_cpu_%a.txt
#SBATCH --array=2

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

NUM_RUNS_PER_SAMPLE=5
sample_index=$((($SLURM_ARRAY_TASK_ID - 1) / $NUM_RUNS_PER_SAMPLE + 1))
run_index=$((($SLURM_ARRAY_TASK_ID - 1) % $NUM_RUNS_PER_SAMPLE + 1))
echo "Executing run $run_index with sample $sample_index..."

repo_dir=$(git rev-parse --show-toplevel)
sample_id_path_1=$repo_dir/code/01_spaceranger/spaceranger_parameters_1v-16v.txt
sample_id_path_2=$repo_dir/code/01_spaceranger/spaceranger_parameters_17v-31v_trimmed.txt

sample_id=$(cat ${sample_id_path_1} ${sample_id_path_2} | awk 'BEGIN {FS=","} {print $1}' | awk "NR==${sample_index}")
in_dir=$repo_dir/processed-data/01_spaceranger/$sample_id/outs
out_dir=$repo_dir/processed-data/17_cell_cell_comm

mkdir -p ${out_dir}

cellnest run \
    --data_name=${sample_id} \
    --model_name=model_${sample_id} \
    --run_id=${run_index} \
    --num_epoch=80000 \
    --model_path=${out_dir}/model/ \
    --embedding_path=${out_dir}/embedding_data/ \
    --training_data=${out_dir}/input_graph/ \
    --manual_seed=yes \
    --seed=$run_index \

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
