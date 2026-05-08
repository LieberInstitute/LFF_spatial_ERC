## Louise Huuki-Myers, May 2026
## Run label transfer to label xenium cells with SingleR

#### Set up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("cell_type_col", "c", "1", "character", "cell type column"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
# opt$cell_type_col <- "cell_type_broad"

data_dir <- here("processed-data", "21_Xenium", "09_xenium_label_transfer_singleR")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("plots", "21_Xenium", "08_xenium_label_transfer_RCTD")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), "- Load data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

spe
# dim: 366 467516 
