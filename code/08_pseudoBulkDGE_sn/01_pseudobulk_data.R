## Louise Huuki-Myers, April 2025
## Pseudobulk single nuc dataset for DGE

library("spatialLIBD")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data")
#if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))
