## Louise Huuki-Myers, April 2025
## Pseudobulk single nuc dataset for DGE

## using pseudobulkDGE_compatable pull request 
# remotes::install_github("LieberInstitute/spatialLIBD", ref = remotes::github_pull(108))
# spatialLIBD          * 1.19.12   2025-04-16 [1] Github (LieberInstitute/spatialLIBD@0e32c65)

library("spatialLIBD")
library("SingleCellExperiment")
library("HDF5Array")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

rowData(sce)

#### run pseudbulk ####
mito_genes <- as.logical(seqnames(sce) == "chrM")
table(mito_genes)

sce_pseudo <- registration_pseudobulk(
    sce,
    var_registration = "cell_type_anno",
    var_sample_id = "sample_id",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL,
    filter_expr = FALSE,
    mito_gene = mito_genes
)

message(Sys.time(), " - Done pseudobulk")
#### Additional edits ####

## drop all NA cols

## add syntactic APOE vars

#### save ####
message(Sys.time(), " - Save")
saveRDS(sce_pseudo, file = here(data_dir, "sce_pseudo_DGE.RDS"))

# slurmjobs::job_single('01_pseudobulk_data', create_shell = TRUE, memory = '100G', command = "Rscript 01_pseudobulk_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
