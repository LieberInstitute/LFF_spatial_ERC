## November 2024, Louise Huuki-Myers
## Compare SVGm vs. Marker clusters

library("spatialLIBD")
library("scater")
library("tidyverse")
library("DeconvoBuddies")
library("bluster")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "18_compare_clusters")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "18_compare_clusters")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### source plotting function ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

#### specify and load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

reducedDimNames(spe)
# [1] "10x_pca"      "10x_tsne"     "10x_umap"     "PCA_p1"       "PCA_p2"       "TSNE"         "UMAP"        
# [8] "HARMONY"      "UMAP.HARMONY" "TSNE.HARMONY"

table(spe$SpD)

offset_clusters <- read_csv(here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_Marker_offset", "BayesSpace_Markers_k09","BayesSpace_Markers_k09", "clusters.csv")) |>
    mutate(key = gsub("(_Br[^_]+)\\1", "\\1", key),
           cluster2 = paste0("Os", cluster))

identical(spe$key, offset_clusters$key)

table(spe$SpD, offset_clusters$cluster2)

#### annotate colors with JACC mat ####
jacc_mat <- linkClustersMatrix(spe$SpD, offset_clusters$cluster2)

jacc_mat_long <- jacc_mat |>
    reshape2::melt() |>
    rename(SpD = Var1, offset = Var2, Jacc = value)
                               

cluster_match <- jacc_mat_long |>
                        group_by(offset) |>
                        arrange(-Jacc) |>
                        slice(1) |>
                        mutate(SpD_match = gsub("~.*?$", "", SpD),
                               SpD_offset = paste0(SpD_match, "~", offset)) |>
                        arrange(SpD_offset)

cluster_match$SpD_offset

#### Match colors ####
## create color pallet

SpD_colors_simple <- metadata(spe)$SpD_colors
names(SpD_colors_simple) <- names(SpD_colors_simple)

SpD_offset_colors <- c(
    "L1~Os4"    = "grey40",   # original L1 blue
    "L1~Os6"    = "#0220DE",   # lighter blue for second L1
    "L2.3~Os1"  = "#FEAF16",
    "L5~Os3"    = "#16FF32",
    "L6~Os5"    = "#178C6D",
    "LD~Os2"    = "#00BCF9",
    "Vasc~Os9"  = "#E05AD2",
    "WM.uf~Os8" = "#E4E1E3",
    "WM~Os7"    = "#581009"
)

offset_clusters <- offset_clusters |> left_join(cluster_match |> select(cluster2 = offset, SpD_offset))

spe$SpD_offset <- offset_clusters$SpD_offset

table(spe$SpD_offset, spe$SpD)

#### Spot plots ####
sample_order <- sort(unique(spe$BrNum))
## plot SVGm with Marker color match

vis_grid_clus(
    spe = spe,
    clustervar = "SpD_offset",
    pdf = here::here(plot_dir, "spe_erc-BayesSpace_offset.pdf"),
    sort_clust = FALSE,
    point_size = 1.2,
    colors = SpD_offset_colors,
    sample_order = sample_order)


#### Jacc heatmaps ####

jacc_mat <- linkClustersMatrix(spe$SpD, spe$SpD_offset)


col_ha <- HeatmapAnnotation(
    df = as.data.frame(list(Markers = colnames(svgVm_jacc_mat$svg11_m11))),
    col = list(
        Markers = SpD_colors$k11_Markers_anno
    ),
    annotation_name_side = "left",
    show_legend = c(FALSE)
)

row_ha <- rowAnnotation(
    df = as.data.frame(list(SVGm = rownames(svgVm_jacc_mat$svg11_m11))),
    col = list(
        SVGm = SpD_colors$k11_SVGm_anno
    ),
    show_legend = c(FALSE)
)

SVGm_count <- rowAnnotation(n_spots = anno_barplot(as.numeric(table(spe$SVGm_k11))))
Mark_count <- columnAnnotation(n_spots = anno_barplot(as.numeric(table(spe$Markers_k11))))

pdf(here(plot_dir, "jacc_mat_k9_offset.pdf"))

print(Heatmap(jacc_mat,
              name = "Correspondence",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              # right_annotation = row_ha,
              # left_annotation = SVGm_count,
              # bottom_annotation = col_ha,
              # top_annotation = Mark_count,
              col = c("black", viridisLite::plasma(100)),
              na_col = "black"
))

dev.off()


#### which cluster has strong correlation with Human Pilot layers? ####

cluster_eval_cor <- do.call("rbind", map2(my_clusters, names(my_clusters), 
                         ~cor_anno[[.x]]$cor_layer$HumanPilot |>
                            reshape2::melt() |>
                             rename(SpD = Var1, HumanPilot = Var2, cor = value) |>
                             mutate(input = .y))) |> 
    left_join(cluster_eval_anno_long)


cluster_eval_cor |> count(anno)


cluster_cor_point <- cluster_eval_cor |>
    # filter(grepl("k11", input)) |>
    ggplot(aes(x = HumanPilot, y = cor, shape = input, color = anno)) +
    geom_point(position = position_jitter(w = 0.2, h = 0), size = 1.5) +
    scale_color_manual(values = c(SpD_colors$k11_Markers_anno, SpD_colors$k11_SVGm_anno, SpD_colors$k09_SVGm_anno)) +
    theme_bw()

ggsave(cluster_cor_point, filename = here(plot_dir, "cluster_cor_point.png"))


#### graph likelihoods together ####

q_tune = map_dfr(list(markers = here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_Marker", "qTune_logliks_Markers.csv"),
         SVGm = here("processed-data", "05_spe_correct_cluster", "05.5_qTune", "qTune_logliks_SVGm.csv")),
    ~read.csv(.x, row.names = 1) |>
        mutate(input = gsub("qTune_logliks_(.*?).csv", "\\1", basename(.x))))

q_tune_line = ggplot(q_tune, 
                     aes(x = q, y = -loglik, color = input)) +
    geom_line() +
    theme_bw()

ggsave(q_tune_line, filename = here(plot_dir, "q_tune_line.png"), height = 5)


