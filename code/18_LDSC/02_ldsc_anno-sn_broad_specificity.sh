#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=02_ldsc_anno-sn_broad_specificity
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-176%20

## Define loops and appropriately subset each variable for the array task ID
all_cluster=(Astro Excit Inhib Macro Micro Oligo OPC Vasc)
cluster=${all_cluster[$(( $SLURM_ARRAY_TASK_ID / 22 % 8 ))]}

all_chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
chr=${all_chr[$(( $SLURM_ARRAY_TASK_ID / 1 % 22 ))]}

## Explicitly pipe script output to a log
log_path=logs/02_ldsc_anno-sn_broad_specificity_${cluster}_${chr}_${SLURM_ARRAY_TASK_ID}.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load LDSC
module load ldsc/1.0.1

## List current modules for reproducibility
module list

BED_DIR="/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/18_LDSC/01_prep_marker_bed/bed_specificity_sn_broad"
OUT_DIR="/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/18_LDSC/LDSC_sn_broad_specificity"

echo "Cluster: "${cluster}
echo "chr: "${chr}

## Annotation
date
echo "Annotation"

python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/make_annot.py \
	--bed-file ${BED_DIR}/${cluster}.bed \
	--bimfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
	--annot-file ${OUT_DIR}/${cluster}/chr.${chr}.annot.gz

## LD score
date
echo "LD score"
python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py \
	 --l2 --thin-annot \
	 --ld-wind-cm 1 \
	 --bfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
	 --anno ${OUT_DIR}/${cluster}/chr.${chr}.annot.gz \
	 --out ${OUT_DIR}/${cluster}/chr.${chr} \
	 --print-snps /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/hapmap3_snps/hm.${chr}.snp


echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

