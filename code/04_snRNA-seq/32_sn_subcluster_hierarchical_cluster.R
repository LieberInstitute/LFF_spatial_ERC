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

# slurmjobs::job_single('32_sn_subcluster_hierarchical_cluster', create_shell = TRUE, memory = '5G', command = "Rscript 32_sn_subcluster_hierarchical_cluster.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


