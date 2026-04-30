## Louise Huuki-Myers, April 2026
## Run label transfer to label xenium cells
## Adapted from https://github.com/LieberInstitute/spatialAmygdala/blob/77670e73360a112945a4b8257557194f9212bdb6/code/Xenium/06_label_transfer/01_RCTD_updated.R
#### Set up ####

library("here")
library("SpatialExperiment")
library("qs2")
library("scater")
library("tidyverse")

## using devel version ‘2.2.1’
library("spacexr")

data_dir <- here("processed-data", "21_Xenium", "08_xenium_label_transfer_RCDT")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "08_xenium_label_transfer_RCDT")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), "- Load data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

spe
# dim: 366 467516 

## Remove any remaining duplicate gene names
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) {
    message("Removing ", sum(dup_genes), " duplicated gene names from spe")
    spe <- spe[!dup_genes, ]
}

any(duplicated(rownames(spe)))


