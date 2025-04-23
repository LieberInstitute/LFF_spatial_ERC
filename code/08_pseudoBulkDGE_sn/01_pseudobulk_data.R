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

# rowData(sce)

## add syntactic APOE vars
sce_pseudo$APOE_syn <- factor(gsub("/", ".", sce_pseudo$APOE))
levels(sce_pseudo$APOE_syn)

sce_pseudo$APOE_carrier_syn <- factor(gsub("\\+", "", sce_pseudo$APOE_carrier))
levels(sce_pseudo$APOE_carrier_syn)

## get mito genes
mito_genes <- as.logical(seqnames(sce) == "chrM")
table(mito_genes)

#### run pseudobulk - cell_type_anno ####
message(Sys.time(), "pseudobulk cell_type_anno")
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

message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))

## Drop cell types w/ not enough pseudobulk samples?
table(sce_pseudo$APOE_carrier, sce_pseudo$cell_type_bulk)

#### Additional edits + Save ####
## drop all NA cols
all_na <- sapply(colData(sce_pseudo), function(x)all(is.na(x)))
colData(sce_pseudo) <- colData(sce_pseudo)[, names(all_na)[!all_na]]

## save
message(Sys.time(), " - Save")
saveRDS(sce_pseudo, file = here(data_dir, "sce_pseudo_DGE.RDS"))

#### run pseudbulk - cell_type_bulk ####
message(Sys.time(), "pseudobulk cell_type_broad")
sce_pseudo <- registration_pseudobulk(
    sce,
    var_registration = "cell_type_bulk",
    var_sample_id = "sample_id",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL,
    filter_expr = FALSE,
    mito_gene = mito_genes
)

message(Sys.time(), " - Done pseudobulk")

message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))

## Drop cell types w/ not enough pseudobulk samples?
table(sce_pseudo$APOE_carrier, sce_pseudo$cell_type_bulk)

#### Additional edits + Save ####
## drop all NA cols
all_na <- sapply(colData(sce_pseudo), function(x)all(is.na(x)))
colData(sce_pseudo) <- colData(sce_pseudo)[, names(all_na)[!all_na]]

## save 
message(Sys.time(), " - Save")
saveRDS(sce_pseudo, file = here(data_dir, "sce_pseudo_DGE_broad.RDS"))


# sce_pseudo <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "sce_pseudo_DGE.RDS"))

# slurmjobs::job_single('01_pseudobulk_data', create_shell = TRUE, memory = '100G', command = "Rscript 01_pseudobulk_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
