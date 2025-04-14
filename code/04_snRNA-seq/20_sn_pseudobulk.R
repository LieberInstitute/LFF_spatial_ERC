## Louise Huuki-Myers, January 2025
## Run spatialLIBD::registration_wrapper to get cell type modeling and pseudobulk data

## Required libraries
library("here")
library("sessioninfo")
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

data_dir <- here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

#### Run Spatial Registration Function ####
cluster_var <- opt$cluster

stopifnot(cluster_var %in% colnames(colData(sce)))

table(sce[[cluster_var]])

message(Sys.time(), " make pseudobulk object")
## I think this needs counts assay
sce_pseudo <- scuttle::aggregateAcrossCells(
    sce,
    DataFrame(
        cluster_var = sce[[cluster_var]]
    )
)

saveRDS(sce_pseudo, file = here(data_dir, sprintf("sce_pseudobulk_only-%s.rds", cluster_var)))

# slurmjobs::job_single('20_sn_pseudobulk', create_shell = TRUE, memory = '100G', command = "Rscript 20_sn_pseudobulk.R -cluster 'cell_type_anno'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
