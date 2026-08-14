## September 2024, Louise Huuki-Myers
## Run spatialLIBD::registration_wrapper to get layer modeling and pseudobulk data


## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("HDF5Array")
library("getopt")

# Import command-line parameters
spec <- matrix(
    c(  "cluster", "i", "append", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(spec)

# cluster_type <- "BayesSpace_SVGm"
cluster_type <- opt$cluster

data_dir <- here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", cluster_type)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## Choose k
k <- as.numeric(
    #   Only one of these environment variables will be defined, so grab the
    #   defined one (handle SGE or SLURM)
    Sys.getenv("SLURM_ARRAY_TASK_ID")
)
stopifnot(!is.na(k))
k_nice <- sprintf("%02d", k)

message("k:", k_nice)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

colnames(colData(spe))

cluster_var <- sprintf("%s_k%02d", cluster_type, k)

## Add syntatic APOE
spe$APOE_syn <- factor(gsub("/", ".", spe$APOE))
levels(spe$APOE_syn)

#### Run Spatial Registration Function ####
message(Sys.time(), " - Running Spatial Registration on: ", cluster_var)
stopifnot(cluster_var %in% colnames(colData(spe)))

modeling_results <-registration_wrapper(
    sce = spe,
    var_registration = cluster_var,
    var_sample_id = "sample_id",
    covars = c("APOE_syn", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    # suffix = k_nice,
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, sprintf("spe_pseudobulk-%s.rds", cluster_var))
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("modeling_results-%s.rds", cluster_var)))

# slurmjobs::job_single('08_model_pseudobulk', create_shell = TRUE, memory = '100G', command = "Rscript 08_model_pseudobulk.R")

# slurmjobs::job_loop('08_model_pseudobulk_Marker', create_shell = TRUE, memory = '100G', loops = list(k = c("2", "9", "15", "20")))

# slurmjobs::array_submit(name = "08_model_pseudobulk", task_ids = c(9), submit = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
        
        
