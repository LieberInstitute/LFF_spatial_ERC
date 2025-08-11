
## Louise Huuki-Myers, Aug 2025
## creat batch jobs for generating annotation files
## based on https://github.com/LieberInstitute/spatialdACC/blob/main/code/17-2_LDSC_enrichment/spatial/ldsc_anno_jobs.pl


#### Set Up ####

library("getopt")
library("here")
library("slurmjobs")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type",
      "mode", "m", "1", "character", "specificity or enrichment"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype <- "Visium"
# opt$datatype <- "sn_broad"
# opt$mode <- "specificity"

print(opt)

## format output dir
out_dir <- here("processed-data", "18_LDSC", sprintf("LDSC_%s_%s", opt$datatype, opt$mode))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_dir

## get bed files/clusters
bed_data_dir <- here("processed-data", "18_LDSC", "01_prep_marker_bed", sprintf("bed_%s_%s", opt$mode, opt$datatype))
bed_files <- list.files(bed_data_dir)
(clusters = gsub(".bed", "", bed_files))

## create sub-output dirs
create_sub_out_dir <- function(name){
    sub_out_dir <- here(out_dir, name)
    if (!dir.exists(sub_out_dir)) dir.create(sub_out_dir, recursive = TRUE)
    
}

purrr::walk(clusters, create_sub_out_dir)
all(clusters %in% list.files(out_dir))


slurmjobs::job_loop(loops = list(cluster = clusters,
                                 chr = as.character(1:22)),
                    create_shell = TRUE,
                    name = sprintf("02_ldsc_anno-%s_%s", opt$datatype, opt$mode),
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

