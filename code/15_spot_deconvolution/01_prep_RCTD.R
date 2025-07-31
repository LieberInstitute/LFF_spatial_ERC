
## Louise Huuki-Myers, July 2025
## Prep data for RCTD spot deconvolution

#### set up ####
library("SpatialExperiment")
library("spacexr")
library("HDF5Array")
library("here")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("cell_type_col", "c", "1", "character", "cell type column"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# plot_dir <- here("plots", "15_spot_deconvolution", "01_prep_RCTD")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "15_spot_deconvolution", "01_prep_RCTD")
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
rctd_data <- createRctd(spe, sce, UMI_min = 2, cell_type_col = opt$cell_type_col)


saveRDS(results_spe, file = here(data_dir, sprintf("rctd_data_%s.rds", opt$cell_type_col)))

slurmjobs::job_loop(loops = list(cell_type_col = c("cell_type_broad", "cell_type_anno")), 
                    create_shell = TRUE, 
                    name = "01_prep_RCTD", 
                    create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
