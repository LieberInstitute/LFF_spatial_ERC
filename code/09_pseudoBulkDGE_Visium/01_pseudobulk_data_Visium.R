## Louise Huuki-Myers, April 2025
## Pseudobulk single nuc dataset for DGE

library("spatialLIBD")
library("SingleCellExperiment")
library("HDF5Array")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

# rowData(spe)

#### run pseudbulk ####
mito_genes <- as.logical(seqnames(spe) == "chrM")
table(mito_genes)

spe_pseudo <- registration_pseudobulk(
    spe,
    var_registration = "vSpD",
    var_sample_id = "sample_id",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL,
    filter_expr = FALSE,
    mito_gene = mito_genes
)

message(Sys.time(), " - Done pseudobulk")

message(sprintf("nrow: %d, ncol: %d", nrow(spe_pseudo), ncol(spe_pseudo)))

## Check SpD frequency
table(spe_pseudo$APOE_carrier, spe_pseudo$vSpD)

#### Additional edits ####
## drop all NA cols
all_na <- sapply(colData(spe_pseudo), function(x)all(is.na(x)))
colData(spe_pseudo) <- colData(spe_pseudo)[, names(all_na)[!all_na]]

## add syntactic APOE vars
spe_pseudo$APOE_syn <- factor(gsub("/", ".", spe_pseudo$APOE))
levels(spe_pseudo$APOE_syn)

spe_pseudo$APOE_carrier_syn <- factor(gsub("\\+", "", spe_pseudo$APOE_carrier))
levels(spe_pseudo$APOE_carrier_syn)

## add E4/E4 variable
spe_pseudo$APOE_E4E4 <- spe_pseudo$APOE == "E4/E4"
table(spe_pseudo$APOE_E4E4)

## add continuous APOE variable
# APOE_num <- data.frame(APOE_num = 0:3)
# rownames(APOE_num) <- c("E2/E2", "E2/E3", "E3/E4", "E4/E4")
# spe_pseudo$APOE_num <- APOE_num[spe_pseudo$APOE, ] 
# table(spe_pseudo$APOE_num, spe_pseudo$APOE)

## add nspots col
spe_pseudo$nspots <- spe_pseudo$ncells

#### Add PCAs ####
message(Sys.time(), " - Run pseudobulk PCA")
reducedDims(spe_pseudo) <- NULL

spe_pseudo <- scater::runPCA(spe_pseudo, 
                         ncomponents = 50,
                         name = "PCA")

#### save ####
message(Sys.time(), " - Save")
saveRDS(spe_pseudo, file = here(data_dir, "spe_pseudo_DGE.RDS"))

# spe_pseudo <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "spe_pseudo_DGE.RDS"))

# slurmjobs::job_single('01_pseudobulk_data_Visium', create_shell = TRUE, memory = '50G', command = "Rscript 01_pseudobulk_data_Visium.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
