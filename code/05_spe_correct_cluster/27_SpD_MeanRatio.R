## Louise Huuki-Myers, April 2025
## Run DeconvoBuddies::get_mean_ratio on snRNA_seq data

## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("DeconvoBuddies")
library("HDF5Array")
library("ggplot2")
library("getopt")
library("dplyr")

# Import command-line parameters
spec <- matrix(
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(spec)
print(opt)

cluster <- opt$cluster

# cluster = "cell_type_fine"

data_dir <- here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "27_SpD_MeanRatio")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
rownames(spe) <- rowData(spe)$gene_name

stopifnot(cluster %in% colnames(colData(spe)))
#### calculate marker stats ####

## Run Mean Ratio
message(Sys.time(), " - Run MeanRatio on ", cluster)
marker_stats_MeanRatio <- get_mean_ratio(
    sce = spe,
    assay_name = "logcounts",
    cellType_col = cluster,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

save(marker_stats_MeanRatio, file = here(data_dir, sprintf("marker_stats_MeanRatio_%s.Rdata", cluster)))


#### plot top markers ####
message(Sys.time(), " - Plots")

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

cluster %in% colnames(colData(spe))

## plot markers
plot_marker_express_ALL(
    spe,
    marker_stats,
    pdf_fn = here(plot_dir, sprintf("SpD_violin_MeanRatio_top10-%s.pdf", cluster)),
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = cluster,
    color_pal = SpD_colors,
    plot_points = FALSE
)


# slurmjobs::job_single('27_SpD_MeanRatio', create_shell = TRUE, memory = '50G', command = "Rscript 27_SpD_MeanRatio.R --cluster 'SpD'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

