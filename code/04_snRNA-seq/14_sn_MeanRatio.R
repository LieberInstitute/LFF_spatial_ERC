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
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

cluster <- opt$cluster

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

## make APOE syntatic
## only ~sample_ID is full rank
# sce$APOE <- gsub("/", "", sce$APOE)
# c("sample_id", "APOE" ,"Sex", "Age","Anc_Afr") %in% colnames(colData(sce))
# mod <- model.matrix(~sample_id + Sex , colData(sce))
# limma::is.fullrank(mod)

message(Sys.time(), " - Run 1vALL on ", cluster)
marker_stats_1vAll <- findMarkers_1vAll(
    sce = sce,
    assay_name = "counts",
    cellType_col = cluster,
    mod = "~sample_id"
)

#### Join and save data ####
message(Sys.time(), " - Done - Join data & save")
## join the two marker_stats tables
marker_stats <- marker_stats_MeanRatio |>
    dplyr::left_join(marker_stats_1vAll, by = join_by(gene, cellType.target))

save(marker_stats, file = here(data_dir, sprintf("MarkerStats_%s.Rdata", cluster)))

# Fix '.' add to cell type names
# marker_stats <- marker_stats |> mutate(cellType.target = gsub("(\\d+)", ".\\1", cellType.target),
#                        cellType.2nd = gsub("(\\d+)", ".\\1", cellType.2nd))

#### plot hockey stick plots & top markers ####
message(Sys.time(), " - Plots")

load(here("processed-data", "04_snRNA-seq", "09_cluster_annotation", "cell_type_colors_allK.Rdata"), verbose = TRUE)

cluster_split <- unlist(strsplit(cluster, split = "_"))
cellType_colors <- cell_type_colors[[cluster_split[[3]]]][[cluster_split[[2]]]]


cluster %in% colnames(colData(sce))

## plot markers
plot_marker_express_ALL(
    sce,
    marker_stats,
    pdf_fn = here(plot_dir, sprintf("sn_violin_MeanRatio_top10-%s.pdf", cluster)),
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = cluster,
    color_pal = cellType_colors,
    plot_points = FALSE
)

## MR vs. 1vALL
hockey_plot <- marker_stats |>
    ggplot(aes(MeanRatio, std.logFC)) +
    geom_point() +
    facet_wrap(~cellType.target) +
    theme_bw() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 1, color = "blue", linetype = "dashed") 

ggsave(hockey_plot, filename = here(plot_dir, sprintf("sn_MRvslogFC-%s.pdf", cluster)), width = 10, height = 10)

## add limit
hockey_plot <- marker_stats |>
    ggplot(aes(MeanRatio, std.logFC)) +
    geom_point() +
    facet_wrap(~cellType.target) +
    theme_bw() +
    xlim(0, 30)+
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 1, color = "blue", linetype = "dashed") 

ggsave(hockey_plot, filename = here(plot_dir, sprintf("sn_MRvslogFC-%s_xlim.pdf", cluster)), width = 10, height = 10)

## free scales
hockey_plot <- marker_stats |>
    ggplot(aes(MeanRatio, std.logFC)) +
    geom_point() +
    facet_wrap(~cellType.target, scales = "free") +
    theme_bw() +
    geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 1, color = "blue", linetype = "dashed") 

ggsave(hockey_plot, filename = here(plot_dir, sprintf("sn_MRvslogFC-%s_free.pdf", cluster)), width = 10, height = 10)

# slurmjobs::job_single('14_sn_MeanRatio', create_shell = TRUE, memory = '100G', command = "Rscript 14_sn_MeanRatio.R -cluster 'ct_broad_k20'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

