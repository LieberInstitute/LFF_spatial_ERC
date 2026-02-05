## Louise Huuki-Myers, Feb 2026
## Cluster Multi-study Oligo

#### Set up ####

library("SingleCellExperiment")
library("HDF5Array")
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")

data_dir <- here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "05_multistudy_Oligo_cluster")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# Import command-line parameters
scec <- matrix(
    c("k", "k", "i", "integer", "k metric for cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

message("Run clustering with k=", opt$k)

#### Load multi-study oligo SCE ####
# message(Sys.time(), " - Load HDF5 sce") ## loads in 3s
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_build", "sce_multistudy_Oligo"))

message(Sys.time(), " - Load Rds sce") ## 
sce <- readRDS(here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_build", "sce_multistudy_Oligo.Rds"))

#### SNN + Walktrap cluster ####
set.seed(2526)

## Build SNN graph
message(Sys.time(), " - running buildSNNGraph: k =", opt$k)
snn.gr <- buildSNNGraph(sce, k = opt$k, use.dimred = "HARMONY")

## Run walk trap clustering
message(Sys.time(), " - running walktrap")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters)

saveRDS(clusters, file = here(data_dir, sprintf("Multistudy_Oligo_cluster_k%i.Rds", opt$k)))

# clusters <- readRDS(file = here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_subcluster", "Multistudy_Oligo_cluster.Rds"))

table(clusters)

# clust.50 <- clusterCells(sce, use.dimred="HARMONY", BLUSPARAM=NNGraphParam(k=50))
# table(clust.50)

message(Sys.time(), " - done walktrap")

#### kmeans clustering #### 
# k = 9
# km <- kmeans(
#     reducedDim(sce, "HARMONY"),
#     centers = 
# )
# table(km$cluster, sce$cell_type_broad)
# clusters <- km$cluster

#### Annotate clusters ####
cluster_anno <- as.data.frame(table(clusters, sce$cell_type_broad)) |>
    pivot_wider(names_from = "Var2", values_from = "Freq") |>
    mutate(oligo_prop = Oligo/(Oligo+OPC),
           cell_type = case_when(oligo_prop > 0.8 ~ "Oligo", 
                                 oligo_prop < 0.2 ~ "OPC",
                                 TRUE ~ "Mix"),
           n = Oligo + OPC) |>
    dplyr::rename(cluster = clusters) |>
    group_by(cell_type) |>
    arrange(cell_type, -n) |>
    mutate(cluster_anno = paste0("multi_", cell_type,".", row_number()))

write.csv(cluster_anno, file = here(data_dir, sprintf("Multistudy_Oligo_cluster_annotation_k%i.csv", opt$k)))

cluster_tab <- data.frame(key = colnames(sce), 
                          cluster = clusters, 
                          cluster_anno = cluster_anno$cluster_anno[match(clusters, cluster_anno$cluster)])

head(cluster_tab)

table(cluster_tab$cluster_anno)

#### save cluster_tab data ####
message(Sys.time() , " - saving data")
write_rds(cluster_tab, file = here(data_dir, sprintf("walktrap_snn_k%02d_subclusters_Multistudy_Oligo.Rds", opt$k)))

## Add to sce data & create color pal
sce$Oligo_anno <- cluster_tab$cluster_anno
other_oligo_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")

sce$metadata$colors <- other_oligo_colors

#### Quality checks ####
pd <- as.data.frame(colData(sce))

colnames(pd)

qc_violin_plot_all <- pd |> 
    select(Oligo_anno, sum, detected, subsets_Mito_percent, doubletScore)  |>
    pivot_longer(!c(Oligo_anno), names_to = "metric") |>
    mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "doubletScore"))) |>
    ggplot() +
    geom_violin(aes(x = Oligo_anno, y = value, fill = Oligo_anno), 
                draw_quantiles = c(0.25, 0.5, 0.75),
                scale = "width") +
    scale_fill_manual(values = other_oligo_colors) +
    theme_bw() +
    facet_grid(metric~., scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(qc_violin_plot_all, filename = here(plot_dir, sprintf("Multistudy_Oligo_k%i_QCmetricViolin.png", opt$k)))


#### Reduced dim plots ####
message(Sys.time(), " - Plot TSNE")
source(here("code", "utils", "my_plot_reduced_dim.R"))

tsne_plot <- my_plot_reduced_dim(sce,
                                 prefix = "other_Oligo",
                                 dimred = "TSNE",
                                 my_var = "Oligo_anno",
                                 var_type = "cat",
                                 save_plot = TRUE,
                                 suffix = sprintf("multistudy_Oligo_k%i", opt$k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE,
                                 color_pal = other_oligo_colors)

# slurmjobs::job_loop(loops = list(k = c("10", "20", "30")),
#                                  name = "05_multistudy_Oligo_cluster",
#                                  create_shell = TRUE,
#                                  create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


