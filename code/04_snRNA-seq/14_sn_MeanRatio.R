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

# Import command-line parameters
scec <- matrix(
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# cluster = "ct_fine_k20"

data_dir <- here("processed-data", "04_snRNA-seq", "14_sn_MeanRatio")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "14_sn_MeanRatio")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

stopifnot(cluster %in% colnames(colData(sce)))

#### calculate marker stats ####

## Run Mean Ratio
message(Sys.time(), " - Run MeanRatio on ", cluster)
marker_stats_MeanRatio <- get_mean_ratio(
    sce = sce,
    assay_name = "logcounts",
    cellType_col = cluster,
    gene_ensembl = "ID",
    gene_name = "Symbol"
)

## Run 1vALL 
message(Sys.time(), " - Run 1vALL on ", cluster)
marker_stats_1vAll <- findMarkers_1vAll(
    sce = sce,
    assay_name = "counts",
    cellType_col = cluster,
    mod = "~BrNum"
)

#### Join and save data ####
message(Sys.time(), " - Done - Join data & save")
## join the two marker_stats tables
marker_stats <- marker_stats_MeanRatio |>
    left_join(marker_stats_1vAll, by = join_by(gene, cellType.target))

save(marker_stats, file = here(data_dir, sprintf("MarkerStats_%s.Rdata", cluster)))


#### plot hockey stick plots & top markers ####

## plot markers
plot_marker_express_ALL(
    sce,
    stats,
    pdf_fn = here(plot_dir, sprintf("sn_violin_MeanRatio_top10-%s.pdf", cluster)),
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = cluster,
    color_pal = NULL,
    plot_points = FALSE
)

## MR vs. 1vALL
hockey_plot <- marker_stats |>
    ggplot(aes(MeanRatio, std.logFC)) +
    geom_point() +
    facet_wrap(~cellType.target) +
    theme_bw()

ggsave(hockey_plot, filename = here(plot_dir, sprintf("sn_MRvslogFC-%s.pdf", cluster)))

# slurmjobs::job_single('14_sn_MeanRatio', create_shell = TRUE, memory = '100G', command = "Rscript 14_sn_MeanRatio.R -cluster 'ct_broad_k20'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

