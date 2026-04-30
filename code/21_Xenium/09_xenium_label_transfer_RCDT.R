## Louise Huuki-Myers, April 2026
## Run label transfer to label xenium cells
## Adapted from https://github.com/LieberInstitute/spatialAmygdala/blob/77670e73360a112945a4b8257557194f9212bdb6/code/Xenium/06_label_transfer/01_RCTD_updated.R
#### Set up ####

library("SpatialExperiment")
library("qs2")
# library("scater")
# library("tidyverse")

# packageVersion("spacexr")
## using devel version ‘2.2.1’
library("spacexr")

library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("cell_type_col", "c", "1", "character", "cell type column"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
opt$cell_type_col <- "cell_type_broad"

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

#### Convert SPE to SpatialRNA (spacexr format) ####

message(Sys.time(), " - Prep RCDT data with refrence: ", opt$cell_type_col)
rctd_data <- createRctd(spe, sce, UMI_min = 2, cell_type_col = opt$cell_type_col, class_df = cell_type_df)


# save
saveRDS(myRCTD, here(processed_dir, "rctd_results_xenium.rds"))