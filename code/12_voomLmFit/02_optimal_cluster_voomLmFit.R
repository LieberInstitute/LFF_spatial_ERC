## Louise Huuki-Myers, July 2025
## Detirmine optimal resolution for DE

#### set-up ####
library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("here")
library("data.table")
library("getopt")
library("tidyverse")
library("sessioninfo")

library("dendextend")
library("dynamicTreeCut")

# Import command-line parameters
scec <- matrix(
    c("celltype", "c", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$celltype <- "Astro"

## set up output dirs
data_dir <- here("processed-data", "12_voomLmFit", "02_optimal_cluster_voomLmFit", opt$celltype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "12_voomLmFit", "02_optimal_cluster_voomLmFit", opt$celltype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

sce <- sce[,sce$cell_type_broad == opt$celltype]

sce$cell_type_anno <- droplevels(sce$cell_type_anno)
table(sce$cell_type_anno)


## hierarchical cluster data
load(here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster", "sn_subcluster_hierarchical_cluster.Rdata"), verbose = TRUE)
# dend
# tree.clusCollapsed
# dist.clusCollapsed

## prune tree to cell types
prune_ct <- labels(dend)[!labels(dend) %in% levels(sce$cell_type_anno)]

dend_ct <- prune(dend, prune_ct)

str(dend_ct)
## get internal node heights
node_heights <- get_nodes_attr(dend_ct, "height")[is.na(get_nodes_attr(dend_ct, "label"))]
cut_heights <- sort(node_heights)[-length(node_heights)]

pdf(here(plot_dir, sprintf("dend_%s.pdf", opt$celltype)))
plot(dend_ct)
abline(h = cut_heights, lty = 2)
dev.off()

hc_ct <- map_dfr(cut_heights, ~cutree(dend_ct, h = .x))
hc_ct$hc.n <- rowMax(as.matrix(hc_ct))

hc_ct_long <- hc_ct |>
    pivot_longer(!hc.n, names_to = "cell_type_anno", values_to = "subgroup") |>
    rowwise() |>
    mutate(hc.ct = paste0("cell_type_hc", hc.n),
           cell_type_hc = paste0(c(opt$celltype, hc.n, subgroup), collapse = ".")
           )

hc_ct_wide <- hc_ct_long |> 
    select(cell_type_anno, cell_type_hc, hc.ct) |>
    pivot_wider(names_from = "hc.ct", values_from = "cell_type_hc")

hc_ct_wide_sce <- hc_ct_wide |>
    column_to_rownames("cell_type_anno")

hc_ct_wide_sce <- hc_ct_wide_sce[sce$cell_type_anno, ]

nrow(hc_ct_wide_sce) == ncol(sce)

colData(sce) <- cbind(colData(sce), hc_ct_wide_sce)

hc_levels <- colnames(hc_ct_wide_sce)

## check tables
map(hc_levels, ~table(sce$cell_type_anno, sce[[.x]]))



