## Louise Huuki-Myers, April 2025
## Run scuttle::aggregateAcrosspells to pseudobulk data

## Required libraries
library("here")
library("sessioninfo")
library("spatialLIBD")
library("HDF5Array")
library("getopt")

# Import command-line parameters
spec <- matrix(
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(spec)
print(opt)

cluster_var <- opt$cluster
# cluster_var <- "SpD"

## prep dir
data_dir <- here("processed-data", "05_spe_correct_cluster", "25_SpD_pseudobulk")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

#### Prep for pseudobulk ####
stopifnot(cluster_var %in% colnames(colData(spe)))

table(spe[[cluster_var]])
class(spe)

## conver to SingleCellExperiment to prevent bug with aggregateAcrossCells
spatialCoords(spe) <- NULL
imgData(spe) <- NULL

spe <- as(spe, "SingleCellExperiment")

## add syntacticly valid version of SpD
if(cluster_var == "SpD"){
    spe$SpD_syn <- factor(gsub("~", "_", spe$SpD), levels = gsub("~", "_", levels(spe$SpD)))
    table(spe$SpD_syn)
    
    cluster_var <- "SpD_syn"
}

#### pseudobulk ####
message(Sys.time(), " make pseudobulk object on ", cluster_var)
class(spe)

spe_pseudo <- scuttle::aggregateAcrossCells(
    spe,
    ids = spe[[cluster_var]]
)

message(sprintf("pseudobulk dims: %d row, %d col", nrow(spe_pseudo), ncol(spe_pseudo)))

## Compute the logcounts
message(Sys.time(), " - normalize expression")
logcounts(spe_pseudo) <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo),
               log = TRUE,
               prior.count = 1
    )

logcounts(spe_pseudo)[1:5,1:5]

if(is(spe_pseudo, "SpatialExperiment")) {
    ## Drop things we don't need
    spatialCoords(spe_pseudo) <- NULL
    imgData(spe_pseudo) <- NULL
}

message(Sys.time(), " - save data")
saveRDS(spe_pseudo, file = here(data_dir, sprintf("spe_pseudobulk_only-%s.rds", cluster_var)))
# spe_pseudo <- readRDS(here(data_dir, sprintf("spe_pseudobulk_only-%s.rds", cluster_var)))

# slurmjobs::job_single('25_SpD_pseudobulk', create_shell = TRUE, memory = '100G', command = "Rscript 25_SpD_pseudobulk.R --cluster 'SpD'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
