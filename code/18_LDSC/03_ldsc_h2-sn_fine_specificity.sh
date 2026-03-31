#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=03_ldsc_h2-sn_fine_specificity
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
log_path=logs/03_ldsc_h2-sn_fine_specificity_${cluster}_${SLURM_ARRAY_TASK_ID}.txt

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

OUT_DIR="/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/18_LDSC/LDSC_sn_fine_specificity"


## H2 
echo "H2"
date
echo "Cluster: "${cluster}
echo "chr: "${chr}

for gwas in adhd alcohol_ldscore Alzheimer_ldscore Alzheimer_ldscore2 Alzheimer_ldscore3 Anorexia_ldscore AUD.EUR_META.NatMed2023 autism_ldscore bmi_ldscore bp_ldscore bp3_ldscore EA_aud_Jun07.txt_ldscore epilepsyAll epilepsyFocal epilepsyGGE GSCAN_AgeSmk_2022_ldscore GSCAN_CigDay_2022_ldscore GSCAN_DrnkWk_2022_ldscore GSCAN_SmkCes_2022_ldscore GSCAN_SmkInit_2022_ldscore height_ldscore insomnia intelligence mdd_ex23andMe_ldscore MDD_ldscore_PGC_UKB_23andme mdd2019edinburgh_ldscore MEGASTROKE.1.AS.EUR_ldscore MEGASTROKE.2.AIS.EUR_ldscore OUD_EA_MVP1_Mar12 OUD_EA_MVP1_MVP2_YP_SAGE_Mar12 OUD_EA_MVP2_Mar12 PASS_Type_2_Diabetes.sumstats PD_ldscore PTSD scz_PGC2_CLOZUK_ldscore scz_PGC3_ldscore smoking_ldscore stroke2022any trait_names_keys UKB_460K.cov_EDU_YEARS.sumstats UKB_460K.mental_NEUROTICISM.sumstats
do
    echo ${gwas}
	
	python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py \
		--h2 /dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/${gwas}.gz \
		--w-ld-chr /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
		--ref-ld-chr ${OUT_DIR}/${cluster}/chr.,/dcs04/lieber/shared/statsgen/LDSC/base/baseline2/baselineLD. \
		--overlap-annot --frqfile-chr /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_frq/1000G.EUR.QC. \
		--out ${OUT_DIR}/${cluster}/${gwas}.out \
		--print-coefficients
done

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

