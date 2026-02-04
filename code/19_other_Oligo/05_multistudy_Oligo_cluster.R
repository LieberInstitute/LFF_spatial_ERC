

#### SNN + Walktrap cluster ####

k = 50

## Build SNN graph
message(Sys.time(), " - running buildSNNGraph: k =", k)
snn.gr <- buildSNNGraph(sce, k = k, use.dimred = "HARMONY")

## Run walk trap clustering
message(Sys.time(), " - running walktrap")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters)

saveRDS(clusters, file = here(data_dir, "Multistudy_Oligo_cluster.Rds"))

# clust.50 <- clusterCells(sce, use.dimred="HARMONY", BLUSPARAM=NNGraphParam(k=50))
# table(clust.50)

message(Sys.time(), " - done walktrap")
cluster_anno <- as.data.frame(table(clusters, sce$cell_type_broad)) |>
    pivot_wider(names_from = "Var2", values_from = "Freq") |>
    mutate(oligo_ratio = Oligo/OPC,
           cell_type = ifelse(oligo_ratio >1, "Oligo", "OPC"),
           n = Oligo + OPC) |>
    dplyr::rename(cluster = clusters) |>
    group_by(cell_type) |>
    arrange(cell_type, -n) |>
    mutate(cluster_anno = paste0(opt$ds_short, "_", cell_type,".", row_number()))

write.csv(cluster_anno, file = here(data_dir, "Multistudy_Oligo_cluster_annotation.csv"))

cluster_tab <- data.frame(key = sce$key, 
                          cluster = clusters, 
                          cluster_anno = cluster_anno$cluster_anno[match(clusters, cluster_anno$cluster)])

head(cluster_tab)

table(cluster_tab$cluster_anno)

## save cluster_tab data
message(Sys.time() , " - saving data")
save(cluster_tab, file = here(data_dir, sprintf("walktrap_snn_k%02d_subclusters_Multistudy_Oligo.Rdata", k)))

## Add to sce data & create color pal
sce$Oligo_anno <- cluster_tab$cluster_anno
other_oligo_colors <- create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")


#### Save data ####

sce$metadata$colors <- other_oligo_colors

saveHDF5SummarizedExperiment(sce, dir = here(data_dir, "sce_multistudy_Oligo"), replace=TRUE)

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

ggsave(qc_violin_plot_all, filename = here(plot_dir, sprintf("other_Oligo_%s_k%i_QCmetricViolin.png", "Multistudy_Oligo", k)))


#### Reduced dim plots ####

message(Sys.time(), " - running TSNE")
sce <- runTSNE(sce, dimred = "HARMONY")

source(here("code", "utils", "my_plot_reduced_dim.R"))

tsne_plot <- my_plot_reduced_dim(sce,
                                 prefix = "other_Oligo",
                                 dimred = "TSNE",
                                 my_var = "Oligo_anno",
                                 var_type = "cat",
                                 save_plot = TRUE,
                                 suffix = sprintf("mulristudy_Oligo_k%i", k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE,
                                 color_pal = other_oligo_colors)
