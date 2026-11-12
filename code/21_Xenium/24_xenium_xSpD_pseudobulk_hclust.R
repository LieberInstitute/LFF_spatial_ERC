## Louise Huuki-Myers, Nov 2026
## Pseudobulk ERC Xenium data by xSpD and run hierarchical clustering
## Adapted from 04_snRNA-seq/20_sn_pseudobulk.R and 04_snRNA-seq/32_sn_subcluster_hierarchical_cluster.R

## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("SingleCellExperiment")
library("qs2")
library("scuttle")
library("edgeR")
library("dendextend")

## Prep directories
data_dir <- here("processed-data", "21_Xenium", "24_xenium_xSpD_pseudobulk_hclust")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "24_xenium_xSpD_pseudobulk_hclust")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

pb_file <- here(data_dir, "spe_xenium_pseudobulk_only-xSpD.rds")

#### Pseudobulk by xSpD (reuse existing output if present) ####
if(file.exists(pb_file)) {
    message(Sys.time(), " - Existing xSpD pseudobulk found, loading: ", basename(pb_file))
    spe_pb <- readRDS(pb_file)
} else {
    message(Sys.time(), " - No xSpD pseudobulk found, running pseudobulk")

    message(Sys.time(), " - Load SpX-level Xenium spe (not filtered to singlets)")
    spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "spe_xenium_bansky.qs2"))
    stopifnot("xSpD" %in% colnames(colData(spe)))

    table(spe$xSpD)
    table(spe$sample_id)

    message(Sys.time(), " - make pseudobulk object")
    ## aggregateAcrossCells keeps per-cell sample_id on a SpatialExperiment, which
    ## cannot map uniquely onto xSpD groups, so aggregate a plain SCE instead
    spe_in <- SingleCellExperiment(
        assays = list(counts = counts(spe)),
        colData = DataFrame(xSpD = spe$xSpD),
        metadata = list(SpX_colors = metadata(spe)$SpX_colors)
    )
    spe_pb <- scuttle::aggregateAcrossCells(
        spe_in,
        DataFrame(
            xSpD = spe_in$xSpD
        )
    )
    
    rm(spe, spe_in)

    colnames(spe_pb) <- spe_pb$xSpD

    message(sprintf("pseudobulk dims: %d row, %d col", nrow(spe_pb), ncol(spe_pb)))

    ## Compute the logcounts from the pseudobulk counts
    message(Sys.time(), " - normalize expression")
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts(spe_pb)))
    logcounts(spe_pb) <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

    ## Drop things we don't need
    if(is(spe_pb, "SpatialExperiment")) {
        spatialCoords(spe_pb) <- NULL
        imgData(spe_pb) <- NULL
    }

    message(Sys.time(), " - save pseudobulk data")
    saveRDS(spe_pb, file = pb_file)
}

#### Perform hierarchical clustering ####
message(Sys.time(), " - Cluster")
dist.clusCollapsed <- dist(t(logcounts(spe_pb)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)

## Save data
message(Sys.time(), " - Save")
save(dend, tree.clusCollapsed, dist.clusCollapsed, file = here(data_dir, "xSpD_hierarchical_cluster.Rdata"))

#### Plot ####
message(Sys.time(), " - Plot")
xspd_colors <- metadata(spe_pb)$SpX_colors[labels(dend)]

pdf(here(plot_dir, "xSpD_dend.pdf"), height = 5)
par(cex = 0.6, font = 2)
dend |>
    set("leaves_col", xspd_colors) |>
    set("leaves_cex", 2) |>
    set("leaves_pch", 19) |>
    plot()
dev.off()

pdf(here(plot_dir, "xSpD_dend_h.pdf"))
par(cex = 0.6, font = 2)
dend |>
    set("leaves_col", xspd_colors) |>
    set("leaves_cex", 2) |>
    set("leaves_pch", 19) |>
    plot(horiz = TRUE)
dev.off()

## Scale point AREA proportionally to number of cells per xSpD
## Largest point has cex = 4
n_cells <- setNames(spe_pb$ncells, colnames(spe_pb))
leaf_n <- n_cells[labels(dend)]
stopifnot(!anyNA(leaf_n), all(leaf_n > 0))
leaf_cex <- 4 * sqrt(leaf_n / max(leaf_n))
legend_counts <- c(1000, 10000, 100000)

pdf(here(plot_dir, "xSpD_dend_ncells.pdf"), height = 5)
par(cex = 0.6, font = 2)
dend |>
    set("leaves_col", xspd_colors) |>
    set("leaves_cex", leaf_cex) |>
    set("leaves_pch", 19) |>
    plot()
legend(
    "topright",
    title = "Number of cells",
    legend = format(legend_counts, big.mark = ",", trim = TRUE),
    pch = 19,
    col = "grey40",
    pt.cex = 4 * sqrt(legend_counts / max(leaf_n)),
    cex = 0.8,
    bty = "n"
)
dev.off()

pdf(here(plot_dir, "xSpD_dend_h_ncells.pdf"))
par(cex = 0.6, font = 2)
dend |>
    set("leaves_col", xspd_colors) |>
    set("leaves_cex", leaf_cex) |>
    set("leaves_pch", 19) |>
    plot(horiz = TRUE)
legend(
    "topleft",
    title = "Number of cells",
    legend = format(legend_counts, big.mark = ",", trim = TRUE),
    pch = 19,
    col = "grey40",
    pt.cex = 4 * sqrt(legend_counts / max(leaf_n)),
    cex = 0.8,
    y.intersp = 1.5,
    bty = "n"
)
dev.off()

# slurmjobs::job_single('24_xenium_xSpD_pseudobulk_hclust', create_shell = TRUE, memory = '50G', command = "Rscript 24_xenium_xSpD_pseudobulk_hclust.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
