#   Convert IF SpatialExperiment, non-IF SpatialExperiment, and snRNA-seq
#   SingleCellExperiment to AnnDatas. Save a copy of the modified SCE as an R
#   object as well. This is a processing step shared by, and prior to, the
#   different deconvolution softwares (tangram, cell2location, SPOTlight).


library("basilisk")
library("SingleCellExperiment")
library("SpatialExperiment")
library("DeconvoBuddies")
library("zellkonverter")
library("sessioninfo")
library("here")
library("HDF5Array")

# cell_group <- "broad" # "broad" or "layer"

#### data prep ####

## SCE
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "sce_objects", "sce_ERC_subcluster"))
sce_in <- here("processed-data", "sce_objects", "sce_ERC_subcluster")

spe_in <- here("processed-data", "spe_objects", "spe_ERC_annotated")

# nonIF_id_path <- here("processed-data", "15_spot_deconvolution", "nonIF_ID_table.csv")
# nonIF_counts_path <- here(
#     "processed-data", "rerun_spaceranger", "{sample_id}", "outs", "spatial",
#     "tissue_spot_counts.csv"
# )

sce_out <- here(
    "processed-data", "15_spot_deconvolution", "01_R_to_python", "sce.h5ad"
)
# sce_r_out <- here(
#     "processed-data", "15_spot_deconvolution", "01_R_to_python",
#     paste0("sce_", cell_group, ".rds")
# )

spe_out <- here(
    "processed-data", "15_spot_deconvolution", "01_R_to_python", "spe.h5ad"
)


marker_object_in <- here(
    "processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio","marker_stats_MeanRatio_cell_type_anno.Rdata"
)

# if (cell_group == "broad") {
#     cell_type_var <- "cellType_broad_hc"
# } else {
#     cell_type_var <- "layer_level"
# }

#  Make sure output directories exist
dir.create(dirname(spe_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(spe_out), recursive = TRUE, showWarnings = FALSE)

####  Functions ####

write_anndata <- function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library("zellkonverter")
                library("reticulate")
                
                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)
                
                #  Write AnnData object to disk
                adata$write(filename = filename)
                
                return()
            },
            env = zellkonverterAnnDataEnv(),
            sce = sce,
            filename = out_path
        )
    )
}

#### Main data conversion ####

#   Load objects
## Load HD5F sce
sce <- HDF5Array::loadHDF5SummarizedExperiment(sce_in)

## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(spe_in)

gc()

####  Check cell counts in spatial object ####
## Use CNmask_dark_blue as nulcei count
spe$count <- spe$CNmask_dark_blue

#   Ensure counts for all spots
if (any(is.na(spe$count))) {
    stop("Did not find cell counts for all non-IF spots.")
}

####   Convert snRNA-seq and spatial R objects to AnnData python objects ####

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
coords <- spatialCoords(spe)
rownames(coords) <- colnames(spe) ## match rownames to colnames(spe)
reducedDims(spe)$spatial <- coords

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) <- rowData(sce)$gene_id

#  convert all objects to Anndatas
message(Sys.time(), " - Converting SCE objects AnnDatas...")
write_anndata(sce, sce_out)

message(Sys.time(), " - Converting SPE objects AnnDatas...")
write_anndata(spe, spe_out)


# #-------------------------------------------------------------------------------
# #   Rank marker genes
# #-------------------------------------------------------------------------------
# 
# #   We won't consider genes that aren't in both the spatial objects
# keep <- (rowData(sce)$gene_id %in% rowData(spe_IF)$gene_id) &
#     (rowData(sce)$gene_id %in% rowData(spe_nonIF)$gene_id)
# sce <- sce[keep, ]
# 
# perc_keep <- 100 * (1 - length(which(keep)) / length(keep))
# print(
#     paste0(
#         "Dropped ", round(perc_keep, 1), "% of potential marker genes ",
#         "that were not present in all the spatial data"
#     )
# )
# 
# print("Running getMeanRatio2 and findMarkers_1vAll to rank genes as markers...")
# marker_stats <- get_mean_ratio2(
#     sce,
#     cellType_col = cell_type_var, assay_name = "logcounts"
# )
# marker_stats_1vall <- findMarkers_1vAll(
#     sce,
#     cellType_col = cell_type_var, assay_name = "logcounts",
#     mod = "~BrNum"
# )
# marker_stats <- left_join(
#     marker_stats, marker_stats_1vall,
#     by = c("gene", "cellType.target")
# )
# 
# saveRDS(marker_stats, marker_object_out)

# slurmjobs::job_single('01_R_to_python', create_shell = TRUE, memory = '25G', command = "Rscript 01_R_to_python.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
