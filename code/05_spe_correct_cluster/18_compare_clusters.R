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
library("readxl")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "18_compare_clusters")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "18_compare_clusters")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### source plotting function ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

#### specify and load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

reducedDimNames(spe)
# [1] "10x_pca"      "10x_tsne"     "10x_umap"     "PCA_p1"       "PCA_p2"       "TSNE"         "UMAP"        
# [8] "HARMONY"      "UMAP.HARMONY" "TSNE.HARMONY"

## create color pallet

SpD_colors <- map(c(9, 10, 11), function(k){
    
    colors <- Polychrome::palette36.colors(k)
    if(k < 3) colors <- colors[1:2]
    names(colors) <- sort(unique(spe[[sprintf("BayesSpace_SVGm_k%02d", k)]]))
    return(colors)
})

my_clusters <- sprintf("k%02d", c(9,10,11))
names(my_clusters) <- my_clusters

names(SpD_colors) <-  my_clusters

#### Add spatial registration data for context ####

# my_clusters <- c(SVGm_k09 = "BayesSpace_SVGm_k09", SVGm_k11 = "BayesSpace_SVGm_k10", SVGm_k11 = "BayesSpace_SVGm_k11")
# load(here("processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC", "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"), verbose = TRUE)
# 
# cluster_eval_anno <- map(my_clusters, 
#                          ~cor_anno[[.x]]$layer_anno$HumanPilot |>
#                              mutate(anno = paste0(layer_label, "~", cluster))
#                              )
# 
# anno_cluster_table <- as.data.frame(map2(my_clusters, cluster_eval_anno, ~.y$anno[match(spe[[.x]], .y$cluster)]))
# head(anno_cluster_table)
# 
# cluster_eval_anno_long <- do.call("rbind", cluster_eval_anno) |>
#     rownames_to_column("input") |>
#     mutate(input = gsub("\\.[0-9]+", "", input)) |>
#     select(input, SpD = cluster, anno)

anno_cluster_table <- map(my_clusters, ~read_xlsx(here("processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC", sprintf("ERC_SpD_spatial_registration_anno_summary_%s.xlsx", .x))))

unique(unlist(map(anno_cluster_table, "anno")))

anno_cluster_table2 <- map2(anno_cluster_table, my_clusters,  ~.x |>
                               mutate(SpD = fct_reorder(paste0(anno, "~", cluster),order)) |>
                               select(!!paste0("BayesSpace_SVGm_", .y) := cluster,
                                      !!paste0("SpD_", .y) := SpD)
                           )

anno_cluster_table3 <- map(anno_cluster_table,  ~.x |>
                               mutate(SpD = fct_reorder(paste0(anno, "~", cluster),order)) |>
                               select(cluster, anno, SpD)
                           )

## Add to colData(spe)
pd <- colData(spe) |> as.data.frame()

pd <- pd[,!grepl("SpD", colnames(pd))]

for(tab in anno_cluster_table2){
    pb <- left_join(pb,tab)
}

colData(spe) <- DataFrame(pb)

table(spe$SpD_k10, spe$BayesSpace_SVGm_k10)

## create color annotation
SpD_colors_V4 <- c("Vasc"      = "#FF56AF",
                   "L1"        = "#47C281",
                   "L2"        = "#41D4EB",
                   # "L3"        = "#0D8278",
                   "L3"        = "#889DF0",
                   "L3_Inhib"  = "#B6686F",
                   "LD"        = "grey80",
                   "L5_LD"     = "#B1C2CE",
                   "L5"        = "#0072CE",
                   "L6"        = "#0A2E5C",
                   "L6a"       = "#24487A",
                   "L6b"       = "#05193B",
                   "WMuf"      = "#F4A460",
                   "WMtz"      = "#E8720C",
                   "WM"        = "#F57A00",
                   "WMd"       = "#581009",
                   "Inhib"     = "#E83E38")

#' Build a color vector for one resolution's annotation table
make_SpD_colors <- function(anno_table, domain_colors = SpD_colors_V4) {
    missing_domains <- setdiff(anno_table$anno, names(domain_colors))
    if (length(missing_domains) > 0) {
        stop("No color defined for domain(s): ", paste(missing_domains, collapse = ", "))
    }
    
    colors_out <- domain_colors[anno_table$anno]
    names(colors_out) <- anno_table$SpD
    colors_out
}

SpD_colors_by_k <- purrr::map(anno_cluster_table3, make_SpD_colors)

#### Compare clusters with Jacquard Matrix ####
svg_marker_pairs <- list(k9_k10 = c("SpD_k09", "SpD_k10"),
                         k9_k11 = c("SpD_k09", "SpD_k11"),
                         k10_k11 = c("SpD_k10", "SpD_k11")
)

svgVm_jacc_mat <- map(svg_marker_pairs, ~linkClustersMatrix(spe[[.x[[1]]]], spe[[.x[[2]]]]))

svgVm_jacc_mat_long <- map(svgVm_jacc_mat, ~.x |> 
                               reshape2::melt() |>
                               rename(SVGm = Var1, Marker = Var2, Jacc = value))

marker_match <- map(svgVm_jacc_mat_long, ~.x |>
                        group_by(SVGm) |>
                        arrange(-Jacc) |>
                        slice(1) |>
                        mutate(Marker_unanno = gsub("~.*?$", "", Marker)) |>
                        arrange(Marker_unanno))

#### Match colors ####

#### Spot plots ####
sample_order <- sort(unique(spe$BrNum))
## plot SVGm with Marker color match

map(my_clusters, ~vis_grid_clus(
    spe = spe,
    clustervar = paste0("SpD_", .x),
    pdf = here::here(plot_dir, sprintf("spe_erc-BayesSpace_SVGm_%s_compare.pdf", .x)),
    sort_clust = FALSE,
    point_size = 1.2,
    colors = SpD_colors_by_k[[.x]],
    sample_order = sample_order)
)

rep_sections_tb <- read.csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "rep_section.csv")) |>
    filter(rep_section)

#### Jacc heatmaps ####
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

pdf(here(plot_dir, "jacc_mat_test.pdf"))

print(Heatmap(svgVm_jacc_mat$svg11_m11,
              name = "Correspondence",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              right_annotation = row_ha,
              # left_annotation = SVGm_count,
              bottom_annotation = col_ha,
              # top_annotation = Mark_count,
              col = c("black", viridisLite::plasma(100)),
              na_col = "black"
))

dev.off()

map2(svgVm_jacc_mat, names(svgVm_jacc_mat), function(jacc.mat, name){
    
    pdf(here(plot_dir, sprintf("jacc_mat_%s.pdf", name)))
    
    print(Heatmap(jacc.mat,
                  name = "Correspondence",
                  # right_annotation = row_ha,
                  # right_annotation = hc_count,
                  # bottom_annotation = col_ha,
                  # bottom_annotation = az_count,
                  col = c("black", viridisLite::plasma(100)),
                  na_col = "black"
    ))
    
    dev.off()
    
})

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


