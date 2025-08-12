#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=02_ldsc_anno-sn_fine_specificity
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-38%20

## Define loops and appropriately subset each variable for the array task ID
all_cluster=(Astro.1 Astro.2 Astro.3 Astro.4 Astro.5 Excit.L2_5.1 Excit.L2_5.2 Excit.L2 Excit.L5_6_NP Excit.L5.1 Excit.L5.2 Excit.L6_CT Excit.L6b Inhib.Chandelier Inhib.Lamp5_Lhx6 Inhib.Pax6 Inhib.Pvalb Inhib.Sst Inhib.Vip Macro Micro.1 Micro.2 Micro.3 Micro.4 Micro.5 Oligo.1 Oligo.2 Oligo.3 Oligo.4 Oligo.5 OPC.1 OPC.2 OPC.3 OPC.4 OPC.5 Vasc.Endo Vasc.PC Vasc.VLMC)
cluster=${all_cluster[$(( $SLURM_ARRAY_TASK_ID / 1 % 38 ))]}

## Explicitly pipe script output to a log
log_path=logs/02_ldsc_anno-sn_fine_specificity_${cluster}_${SLURM_ARRAY_TASK_ID}.txt

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

BED_DIR="/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/18_LDSC/01_prep_marker_bed/bed_specificity_sn_fine"
OUT_DIR="/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/18_LDSC/LDSC_sn_fine_specificity"

echo "Cluster: "${cluster}

for chr in $(seq 1 22);
do
    echo "Chr: " ${chr}
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
	
done

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

