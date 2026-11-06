## Louise Huuki-Myers, May 2025
## Run hierarchical clustering on sn subclustering to examine different resolutions

library("SingleCellExperiment")
library("dendextend")
library("dynamicTreeCut")
library("ggraph")
library("here")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### load pseudobulk data ####
sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk", "sce_ERC_subcluster_pseudobulk_only-cell_type_anno.rds"))
colnames(sce_pb)


#### Perform hierarchical clustering ####
message(Sys.time(), " - Cluster Again")
dist.clusCollapsed <- dist(t(logcounts(sce_pb)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

# dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)

## Rotate oligodendrocyte branches into the desired display order
oligo_order <- c("Oligo.3", "Oligo.5", "Oligo.4", "Oligo.1", "Oligo.2")

# Keep the other cell types in their existing positions
leaf_order <- tree.clusCollapsed$labels[tree.clusCollapsed$order]
stopifnot(all(oligo_order %in% leaf_order))
leaf_order[leaf_order %in% oligo_order] <- oligo_order

# Update the hclust object, then regenerate the plotting dendrogram
tree.clusCollapsed <- dendextend::rotate(
    tree.clusCollapsed,
    order = leaf_order
)
dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)


## Save data
message(Sys.time(), " - Save")
save(dend, tree.clusCollapsed, dist.clusCollapsed, file = here(data_dir, "sn_subcluster_hierarchical_cluster.Rdata"))
# load(here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster", "sn_subcluster_hierarchical_cluster.Rdata"))


#### plot ####
message(Sys.time(), " - Plot")
load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)


broad_clus = cutree(dend, 7)

# Print for future reference
pdf(here(plot_dir, "sn_subcluster_dend.pdf"), height = 5)
par(cex = 0.6, font = 2)

## cell type colors on leaves & branches
dend |> 
    set("leaves_col", cell_type_colors$anno[labels(dend)]) |> 
    set("leaves_cex", 2) |> 
    set("leaves_pch", 19) |> 
    set("branches_k_color", k = 7,
        value = cell_type_colors$broad[c("Inhib", "Excit", "Micro", "Vasc", "Oligo", "Astro", "OPC")]) |>
    plot()

# abline(v = 500, lty = 2)
dev.off()


pdf(here(plot_dir, "sn_subcluster_dend_h.pdf"))
par(cex = 0.6, font = 2)

dend |> 
    set("leaves_col", cell_type_colors$anno[labels(dend)]) |> 
    set("leaves_cex", 2) |> 
    set("leaves_pch", 19) |> 
    set("branches_k_color", k = 7,
        value = cell_type_colors$broad[c("Inhib", "Excit", "Micro", "Endo", "Oligo", "Astro", "OPC")]) |>
    plot(horiz = TRUE)

abline(v = 500, lty = 2)
dev.off()

## Match nucleus counts to dendrogram leaf order
n_cells <- setNames(sce_pb$ncells, colnames(sce_pb))
leaf_n <- n_cells[labels(dend)]

stopifnot(
    length(leaf_n) == length(labels(dend)),
    !anyNA(leaf_n),
    all(leaf_n > 0)
)

## Scale point AREA proportionally to nucleus count
## Largest point has cex = 4
leaf_cex <- 4 * sqrt(leaf_n / max(leaf_n))

legend_counts <- c(100, 1000, 10000)

pdf(here(plot_dir, "sn_subcluster_dend_ncells.pdf"), height = 5)
par(cex = 0.6, font = 2)

## cell type colors on leaves & branches
dend |> 
    set("leaves_col", cell_type_colors$anno[labels(dend)]) |> 
    set("leaves_cex", leaf_cex) |>
    set("leaves_pch", 19) |> 
    set("branches_k_color", k = 7,
        value = cell_type_colors$broad[c("Inhib", "Excit", "Micro", "Vasc", "Oligo", "Astro", "OPC")]) |>
    plot()

legend(
    "topright",
    title = "Number of nuclei",
    legend = format(legend_counts, big.mark = ",", trim = TRUE),
    pch = 19,
    col = "grey40",
    pt.cex = 4 * sqrt(legend_counts / max(leaf_n)),
    cex = 0.8,
    bty = "n"
)

# abline(v = 500, lty = 2)
dev.off()


pdf(here(plot_dir, "sn_subcluster_dend_h_ncells.pdf"))
par(cex = 0.6, font = 2)

dend |> 
    set("leaves_col", cell_type_colors$anno[labels(dend)]) |> 
    set("leaves_cex", leaf_cex) |>
    set("leaves_pch", 19) |> 
    set("branches_k_color", k = 7,
        value = cell_type_colors$broad[c("Inhib", "Excit", "Micro", "Endo", "Oligo", "Astro", "OPC")]) |>
    plot(horiz = TRUE)

abline(v = 500, lty = 2)

legend(
    "topleft",
    title = "Number of nuclei",
    legend = format(legend_counts, big.mark = ",", trim = TRUE),
    pch = 19,
    col = "grey40",
    pt.cex = 4 * sqrt(legend_counts / max(leaf_n)),
    cex = 0.8,
    y.intersp = 1.5,
    bty = "n"
)

dev.off()


# slurmjobs::job_single('32_sn_subcluster_hierarchical_cluster', create_shell = TRUE, memory = '5G', command = "Rscript 32_sn_subcluster_hierarchical_cluster.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


