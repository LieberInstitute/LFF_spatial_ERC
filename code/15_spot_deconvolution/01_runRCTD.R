
## Louise Huuki-Myers, July 2025
## Run RCTD spot deconvolution

#### set up ####
library("ggplot2")
library("SpatialExperiment")
library("SummarizedExperiment")
library("spacexr")
library("HDF5Array")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "15_spot_deconvolution", "01_runRCTD")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "15_spot_deconvolution", "01_runRCTD")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load reference data ####
## Load HD5F spe
message(Sys.time(), " - Load sn ref data")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "sce_objects", "sce_ERC_subcluster"))

# Row names should be unique gene names
rowData(sce)
stopifnot(!any(duplicated(rownames(sce))))

# column names should be unique cell barcodes
colnames(sce) <- sce$Barcode

# colData(sce)
stopifnot(!any(duplicated(colnames(sce))))

# copy sum umi col to nUMI
sce$nUMI <- sce$sum

## cell type
table(sce$cell_type_anno)

#### load spatial query data ####
message(Sys.time(), " - Load Spatial data")
## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))

head(spatialCoords(spe))

# copy sum umi col to nUMI
spe$nUMI <- spe$sum_umi

# Row names should be unique gene names
rowData(spe)
stopifnot(!any(duplicated(rownames(spe))))

# column names should be unique spot barcodes
# colData(spe)
stopifnot(!any(duplicated(colnames(spe))))

# Examine pixels. (optional)
print(dim(assay(spe)))

#### prep RTCD data ####
message(Sys.time(), " - Prep RCDT data")
rctd_data <- createRctd(spe, sce, UMI_min = 2, cell_type_col = "cell_type_anno")

#### Run RTCD  ####
message(Sys.time(), " - Run RCDT")
## max_multi_types = median(spe$CNmask_dark_blue)
results_spe <- runRctd(rctd_data, rctd_mode = "multi", max_multi_types = 5)

message(Sys.time(), " - Done save")
saveRDS(results_spe, file = here(data_dir, "RCDT_est_prop.rds"))

# slurmjobs::job_single('01_runRCTD', create_shell = TRUE, memory = '200G', command = "Rscript 01_runRCTD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
