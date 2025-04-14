## Louise Huuki-Myers, April 2025
## Run scuttle::aggregateAcrossCells to pseudobulk data

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

colnames(sce_pseudo) <- sce_pseudo$cluster_var

## Drop lowly-expressed genes
# message(Sys.time(), " drop lowly expressed genes")
# keep_expr <-
#     edgeR::filterByExpr(sce_pseudo, group = sce_pseudo$cluster_var)
# sce_pseudo <- sce_pseudo[which(keep_expr), ]

message(sprintf("pseudobulk dims: %d row, %d col", nrow(sce_pseudo), ncol(sce_pseudo)))

## Compute the logcounts
message(Sys.time(), " - normalize expression")
logcounts(sce_pseudo) <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
               log = TRUE,
               prior.count = 1
    )

logcounts(sce_pseudo)[1:5,1:5]

if(is(sce_pseudo, "SpatialExperiment")) {
    ## Drop things we don't need
    spatialCoords(sce_pseudo) <- NULL
    imgData(sce_pseudo) <- NULL
}

message(Sys.time(), " - save data")
saveRDS(sce_pseudo, file = here(data_dir, sprintf("sce_pseudobulk_only-%s.rds", cluster_var)))
sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_only-%s.rds", cluster_var)))

# slurmjobs::job_single('20_sn_pseudobulk', create_shell = TRUE, memory = '100G', command = "Rscript 20_sn_pseudobulk.R -cluster 'cell_type_anno'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
