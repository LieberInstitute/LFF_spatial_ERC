## Louise Huuki-Myers, January 2025
## Run DeconvoBuddies::get_mean_ratio on snRNA_seq data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("DeconvoBuddies")
library("HDF5Array")
library("ggplot2")
library("getopt")
library("dplyr")

# Import command-line parameters
scec <- matrix(
    c("cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

cluster <- opt$cluster

# cluster <- "cell_type_anno"

data_dir <- here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "34_sn_subcluster_MeanRatio")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

stopifnot(cluster %in% colnames(colData(sce)))
#### calculate marker stats ####

## Run Mean Ratio
message(Sys.time(), " - Run MeanRatio on ", cluster)
marker_stats_MeanRatio <- get_mean_ratio(
    sce = sce,
    assay_name = "logcounts",
    cellType_col = cluster,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

save(marker_stats_MeanRatio, file = here(data_dir, sprintf("marker_stats_MeanRatio_%s.Rdata", cluster)))
# load(here(data_dir, sprintf("marker_stats_MeanRatio_%s.Rdata", cluster)))

#### export top MR genes ####
marker_stats_MeanRatio |> 
    arrange(cellType.target) |>
    filter(MeanRatio.rank <= 50) |>
    write.csv(here(data_dir, sprintf("sn_top_MeanRatio_%s.csv", cluster)), row.names = FALSE)


#### plot hockey stick plots & top markers ####
message(Sys.time(), " - Plots")

cell_type_colors <- metadata(sce)$cell_type_colors

cluster %in% colnames(colData(sce))

## plot markers
plot_marker_express_ALL(
    sce,
    marker_stats_MeanRatio,
    pdf_fn = here(plot_dir, sprintf("sn_violin_MeanRatio_top10-%s.pdf", cluster)),
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = cluster,
    color_pal = cell_type_colors$anno,
    plot_points = FALSE
)


# slurmjobs::job_single('34_sn_subcluster_MeanRatio', create_shell = TRUE, memory = '100G', command = "Rscript 34_sn_subcluster_MeanRatio.R -cluster 'cell_type_anno'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

