
## Louise Huuki-Myers, Aug 2025
## creat batch jobs for generating annotation files
## based on https://github.com/LieberInstitute/spatialdACC/blob/main/code/17-2_LDSC_enrichment/spatial/ldsc_anno_jobs.pl


#### Set Up ####

library("getopt")
library("here")
library("slurmjobs")

# Import command-line parameters
# scec <- matrix(
#     c("datatype", "d", "1", "character", "Data type",
#       "mode", "m", "1", "character", "specificity or enrichment"),
#     ncol = 5, byrow = TRUE
# )
# opt <- getopt(scec)
# print(opt)

## test
# datatype <- "Visium"
# datatype <- "sn_fine"
# datatype <- "sn_broad"

mode <- "specificity"

## format output dir
out_dir <- here("processed-data", "18_LDSC", sprintf("LDSC_%s_%s", datatype, mode))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_dir

#### 2. LDSC annotation ####
## get bed files/clusters
bed_data_dir <- here("processed-data", "18_LDSC", "01_prep_marker_bed", sprintf("bed_%s_%s", mode, datatype))
bed_files <- list.files(bed_data_dir)
(clusters = gsub(".bed", "", bed_files))

## create sub-output dirs
create_sub_out_dir <- function(name){
    sub_out_dir <- here(out_dir, name)
    if (!dir.exists(sub_out_dir)) dir.create(sub_out_dir, recursive = TRUE)
    }

purrr::walk(clusters, create_sub_out_dir)
all(clusters %in% list.files(out_dir))

# for i in $(seq 1 22);
# do
# echo $i
# done

slurmjobs::job_loop(loops = list(cluster = clusters),
                    create_shell = TRUE,
                    name = sprintf("02_ldsc_anno-%s_%s", datatype, mode),
                    create_script = FALSE)

## check reference inputs on JHPCE
make_annot="/dcs04/lieber/shared/statsgen/LDSC/base/scripts/make_annot.py"
file.exists(make_annot)

referenceDir="/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/"
dir.exists(referenceDir)

# python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/make_annot.py 
# --bed-file bedfiles/L1.bed --bimfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bim --annot-file ./out_L1/chr.1.annot.gz

# cat(sprintf("python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/make_annot.py \
#         --bed-file %s/Astro.bed \
#         --bimfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bim \
#         --annot-file %s/Astro/chr.1.annot.gz",
# )

cat(
    sprintf(
'## Load LDSC
module load ldsc/1.0.1

## List current modules for reproducibility
module list

BED_DIR="%s"
OUT_DIR="%s"

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
	
done',
bed_data_dir,
out_dir
))


#### LD Score ####

#### LDSC H2 ####

ldsc = "/dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py"
baseline2 = "/dcs04/lieber/shared/statsgen/LDSC/base/baseline2/baselineLD."
weights = "/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
freq = "/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_frq/1000G.EUR.QC."

gwas_datasets <- gsub(".gz", "", list.files("/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/"))
##loop over gwas datasets in sh file

slurmjobs::job_loop(loops = list(cluster = clusters),
                    create_shell = TRUE,
                    name = sprintf("03_ldsc_h2-%s_%s", datatype, mode),
                    create_script = FALSE)

cat(sprintf(
    '## Load LDSC
module load ldsc/1.0.1

## List current modules for reproducibility
module list

OUT_DIR="%s"


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
',
    out_dir
))


