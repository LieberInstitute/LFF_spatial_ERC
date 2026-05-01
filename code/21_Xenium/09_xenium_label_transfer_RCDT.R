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
# opt$cell_type_col <- "cell_type_broad"

data_dir <- here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCDT")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("plots", "21_Xenium", "08_xenium_label_transfer_RCDT")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), "- Load data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

spe
# dim: 366 467516 

## create cell type class table for fine cell type
cell_type_df <- NULL
if(opt$cell_type_col == "cell_type_anno"){
    
    message("Use class_df")
    cell_type_levels <- levels(sce$cell_type_anno)
    
    cell_type_df <- data.frame(class = jaffelab::ss(cell_type_levels, "\\."))
    rownames(cell_type_df) <- cell_type_levels
    
}

#### Convert SPE to SpatialRNA (spacexr format) ####
message(Sys.time(), " - Convert SPE to Reference")

## Remove any remaining duplicate gene names
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) {
    message("Removing ", sum(dup_genes), " duplicated gene names from spe")
    spe <- spe[!dup_genes, ]
}

## Extract raw counts
spatial_counts <- counts(spe)

## Final safety check on the extracted matrix
stopifnot("Duplicate rownames still present in spatial_counts" =
              !any(duplicated(rownames(spatial_counts))))

## Extract spatial coordinates as a data.frame with columns x, y
coords <- as.data.frame(spatialCoords(spe))
if (ncol(coords) == 2) {
    colnames(coords) <- c("x", "y")
} else {
    stop("Expected 2 coordinate columns in spatialCoords(spe)")
}

## Ensure rownames on coords match colnames of counts
barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

## nUMI per spot/cell
nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- barcodes

## Sanity check
stopifnot(
    "Barcode mismatch" =
        identical(rownames(coords), colnames(spatial_counts)) &&
        identical(names(nUMI_spatial), colnames(spatial_counts))
)

## Create SpatialRNA object
puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# rownames(spe)

#### Convert SCE to Reference (spacexr format) ####
message(Sys.time(), " - Convert SCE to Reference with cell_types = ", opt$cell_type_col)
rownames(sce) <- rowData(sce)$gene_name

colnames(sce) <- sce$Barcode

## Remove duplicate genes from reference
dup_ref <- duplicated(rownames(sce))
if (any(dup_ref)) {
    message("Removing ", sum(dup_ref), " duplicated gene names from sce")
    sce <- sce[!dup_ref, ]
}

ref_counts <- counts(sce)
ref_counts[1:5, 1:5]


## Cell type labels
cell_types <- sce[[opt$cell_type_col]]
names(cell_types) <- colnames(sce)

## nUMI per cell
nUMI_ref <- Matrix::colSums(ref_counts)
names(nUMI_ref) <- colnames(sce)

## Create Reference object
reference <- Reference(ref_counts, cell_types, nUMI_ref)

#### create object & run RCTD ####
message(Sys.time(), " - create RCDT data")

rctd_data <- create.RCTD(spatialRNA = spe, 
                         reference = sce, 
                         max_cores = 4,
                         class_df = cell_type_df)


message(Sys.time(), "Run RCDT")
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet", cell_type_df)

message(Sys.time(), " - Save qs2")
qs2::qs_save(myRCTD, here(data_dir, "rctd_results_xenium.qs2"))

# slurmjobs::job_single('09_xenium_label_transfer_RCDT', create_shell = TRUE, memory = '200G', command = "Rscript 09_xenium_label_transfer_RCDT.R --cell_type_col cell_type_anno")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()



