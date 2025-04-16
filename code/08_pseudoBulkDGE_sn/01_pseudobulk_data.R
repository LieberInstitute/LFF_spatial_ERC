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

message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))

## Drop cell types w/ not enough pseudobulk samples?
table(sce_pseudo$APOE_carrier, sce_pseudo$cell_type_anno)
# Astro.1 Astro.2 Endo Micro.1 Micro.2 Oligo.1 Oligo.2 OPC Excit.L2 Excit.L2_5.1 Excit.L2_5.2 Excit.L5.1 Excit.L5.2
# E2+      14      14   14      14      14      14      14  14       10           14            6         11          4
# E4+      17      17   16      17      17      17      17  17       15           17            9         14          6
# 
# Excit.L5_6_NP Excit.L6_CT Excit.L6b Inhib.Pax6 Inhib.Lamp5_Lhx6 Inhib.Pvalb Inhib.Vip Inhib.Chandelier Inhib.Sst
# E2+             1           7        10          0               12          11        14                4        12
# E4+             5          10         9          2               15          16        17               11        15

#### Additional edits ####
## drop all NA cols
all_na <- sapply(colData(sce_pseudo), function(x)all(is.na(x)))
colData(sce_pseudo) <- colData(sce_pseudo)[, names(all_na)[!all_na]]

## add syntactic APOE vars
sce_pseudo$APOE_syn <- factor(gsub("/", ".", sce_pseudo$APOE))
levels(sce_pseudo$APOE_syn)

sce_pseudo$APOE_carrier_syn <- factor(gsub("\\+", "", sce_pseudo$APOE_carrier))
levels(sce_pseudo$APOE_carrier_syn)

#### save ####
message(Sys.time(), " - Save")
saveRDS(sce_pseudo, file = here(data_dir, "sce_pseudo_DGE.RDS"))

# sce_pseudo <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "sce_pseudo_DGE.RDS"))

# slurmjobs::job_single('01_pseudobulk_data', create_shell = TRUE, memory = '100G', command = "Rscript 01_pseudobulk_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
