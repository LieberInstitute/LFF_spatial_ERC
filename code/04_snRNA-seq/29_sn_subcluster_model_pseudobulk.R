## Louise Huuki-Myers, MAy 2025
## Run spatialLIBD::registration_wrapper to get cell type subcluster modeling and pseudobulk data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("spatialLIBD")
library("HDF5Array")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

data_dir <- here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

#### Run Spatial Registration Function ####
cluster_var <- opt$cluster

## make APOE syntatic
colData(sce)[,c(cluster_var, "sample_id", "APOE", "Sex", "Age", "Anc_Afr")]
sce$APOE <- gsub("/", "", sce$APOE)

message(Sys.time(), " - Running Spatial Registration on: ", cluster_var)
stopifnot(cluster_var %in% colnames(colData(sce)))

table(sce[[cluster_var]], sce$sample_id)

modeling_results <-registration_wrapper(
    sce = sce,
    var_registration = cluster_var,
    var_sample_id = "sample_id",
    covars = c("APOE", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, sprintf("sce_subcluster_pseudobulk-%s.rds", cluster_var))
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("sce_subcluster_modeling_results-%s.rds", cluster_var)))

# slurmjobs::job_single('29_sn_subcluster_model_pseudobulk', create_shell = TRUE, memory = '100G', command = "Rscript 29_sn_subcluster_model_pseudobulk.R --cluster 'cell_type_anno'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
        
        
