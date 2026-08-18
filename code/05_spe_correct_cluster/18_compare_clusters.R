## November 2024, Louise Huuki-Myers
## Compare SVGm clusters across resolutions (k=9, k=10, k=11)

library("spatialLIBD")
library("tidyverse")
library("bluster")
library("ComplexHeatmap")
library("patchwork")
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

my_clusters <- sprintf("k%02d", c(9, 10, 11))
names(my_clusters) <- my_clusters

#### Load spatial registration annotations ####
anno_cluster_table <- map(
    my_clusters,
    ~ read_xlsx(here(
        "processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC",
        sprintf("ERC_SpD_spatial_registration_anno_summary_%s.xlsx", .x)
    ))
)

unique(unlist(map(anno_cluster_table, "anno")))

anno_cluster_table2 <- map2(
    anno_cluster_table, my_clusters,
    ~ .x |>
        mutate(SpD = fct_reorder(paste0(anno, "~", cluster), order)) |>
        select(
            !!paste0("BayesSpace_SVGm_", .y) := cluster,
            !!paste0("SpD_", .y) := SpD
        )
)

anno_cluster_table3 <- map(
    anno_cluster_table,
    ~ .x |>
        mutate(SpD = fct_reorder(paste0(anno, "~", cluster), order)) |>
        select(cluster, anno, SpD)
)

## Add SpD_k09 / SpD_k10 / SpD_k11 to colData(spe)
pd <- colData(spe) |> as.data.frame()
pd <- pd[, !grepl("SpD", colnames(pd))]

for (tab in anno_cluster_table2) {
    pd <- left_join(pd, tab)
}

colData(spe) <- DataFrame(pd)

table(spe$SpD_k10, spe$BayesSpace_SVGm_k10)

#### Build color palette, keyed by domain name ####
SpD_colors_V4 <- c(
    "Vasc"      = "#FF56AF",
    "L1"        = "#47C281",
    "L2"        = "#41D4EB",
    # "L3"      = "#0D8278",
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
    "Inhib"     = "#E83E38"
)

#' Build a color vector for one resolution's annotation table, keyed by SpD
make_SpD_colors <- function(anno_table, domain_colors = SpD_colors_V4) {
    missing_domains <- setdiff(anno_table$anno, names(domain_colors))
    if (length(missing_domains) > 0) {
        stop("No color defined for domain(s): ", paste(missing_domains, collapse = ", "))
    }

    colors_out <- domain_colors[anno_table$anno]
    names(colors_out) <- anno_table$SpD
    colors_out
}

## one color vector per resolution: SpD_colors_by_k$k09, $k10, $k11
SpD_colors_by_k <- purrr::map(anno_cluster_table3, make_SpD_colors)

#### Spot plots, colored by SpD ####
sample_order <- sort(unique(spe$BrNum))

map(my_clusters, ~ vis_grid_clus(
    spe = spe,
    clustervar = paste0("SpD_", .x),
    pdf = here(plot_dir, sprintf("spe_erc-BayesSpace_SVGm_%s_compare.pdf", .x)),
    sort_clust = FALSE,
    point_size = 1.2,
    colors = SpD_colors_by_k[[.x]],
    sample_order = sample_order
))

#### Representative section grid, for main figures ####
## 4 representative sections (rep_section_4) as columns x k09/k10/k11 as rows.
## Legend is only shown on the right-most sample of each row, since colors
## are consistent within a row (SpD_colors_by_k[[k]]).

rep_sections_tb <- read.csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "rep_section.csv"))

rep_samples_4 <- rep_sections_tb |>
    filter(rep_section_4) |>
    pull(sample_id)

cluster_plots_by_k <- map(my_clusters, function(k) {
    cluster_row_plots <- imap(rep_samples_4, function(s, i) {
        is_last <- i == length(rep_samples_4)

        vis_clus_plot <- vis_clus(
            spe = spe,
            point_size = 1.5,
            colors = SpD_colors_by_k[[k]],
            sampleid = s,
            clustervar = paste0("SpD_", k)
        ) +
            labs(title = s, y = if (i == 1) k else NULL) +
            theme(
                legend.position = if (is_last) "right" else "none", ## only right-most sample keeps its legend
                axis.title.x = element_blank(),
                text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5)
            )

        vis_clus_plot
    })

    Reduce("+", cluster_row_plots) + plot_layout(nrow = 1)
})

cluster_grid <- Reduce("/", cluster_plots_by_k)

ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.pdf"), width = 18, height = 9)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.png"), width = 18, height = 9)

#### Compare resolutions with a Jaccard-style correspondence matrix ####
k_pairs <- list(
    k09_k10 = c("SpD_k09", "SpD_k10"),
    k09_k11 = c("SpD_k09", "SpD_k11"),
    k10_k11 = c("SpD_k10", "SpD_k11")
)

jacc_mat_list <- map(k_pairs, ~ linkClustersMatrix(spe[[.x[[1]]]], spe[[.x[[2]]]]))

jacc_mat_long <- map(jacc_mat_list, ~ .x |>
    reshape2::melt() |>
    rename(SpD_row = Var1, SpD_col = Var2, Jacc = value))

## for each cluster in the first resolution, its single best-matching cluster
## in the second resolution
best_match <- imap_dfr(jacc_mat_long, ~ .x |>
    group_by(SpD_row) |>
    arrange(-Jacc) |>
    slice(1) |>
    ungroup() |>
    mutate(
        match_domain = gsub("~.*?$", "", SpD_col),
        pair = .y
    ) |>
    arrange(match_domain))

write.csv(best_match, file = here(data_dir, "SpD_best_jaccard_match.csv"), row.names = FALSE)

#### Jaccard heatmaps, annotated with each resolution's own SpD colors ####

#' Build a row (or column) annotation for a Jaccard heatmap, named after the
#' resolution it represents and colored from that resolution's SpD_colors_by_k
#' entry, so k09/k10/k11 each keep their own consistent domain colors.
make_jacc_annotation <- function(labels, k_label, colors_by_k, side = c("row", "column")) {
    side <- match.arg(side)
    ann_df <- setNames(data.frame(labels, stringsAsFactors = FALSE), k_label)
    ann_col <- setNames(list(colors_by_k[[k_label]]), k_label)

    if (side == "row") {
        rowAnnotation(df = ann_df, col = ann_col, show_legend = FALSE, annotation_name_side = "top")
    } else {
        HeatmapAnnotation(df = ann_df, col = ann_col, show_legend = FALSE)
    }
}

plot_jacc_heatmap <- function(jacc_mat, row_k, col_k, colors_by_k) {
    row_ha <- make_jacc_annotation(rownames(jacc_mat), row_k, colors_by_k, side = "row")
    col_ha <- make_jacc_annotation(colnames(jacc_mat), col_k, colors_by_k, side = "column")

    Heatmap(
        jacc_mat,
        name = "Correspondence",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = row_ha,
        bottom_annotation = col_ha,
        col = c("black", viridisLite::plasma(100)),
        na_col = "black",
        row_title = row_k,
        column_title = col_k
    )
}

iwalk(jacc_mat_list, function(jacc_mat, pair_name) {
    row_k <- gsub("SpD_", "", k_pairs[[pair_name]][[1]])
    col_k <- gsub("SpD_", "", k_pairs[[pair_name]][[2]])

    pdf(here(plot_dir, sprintf("jacc_mat_%s.pdf", pair_name)))
    print(plot_jacc_heatmap(jacc_mat, row_k, col_k, SpD_colors_by_k))
    dev.off()
})

#### graph qTune likelihoods together ####
q_tune <- map_dfr(
    list(
        markers = here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_Marker", "qTune_logliks_Markers.csv"),
        SVGm = here("processed-data", "05_spe_correct_cluster", "05.5_qTune", "qTune_logliks_SVGm.csv")
    ),
    ~ read.csv(.x, row.names = 1) |>
        mutate(input = gsub("qTune_logliks_(.*?).csv", "\\1", basename(.x)))
)

q_tune_line <- ggplot(q_tune, aes(x = q, y = -loglik, color = input)) +
    geom_line() +
    theme_bw()

ggsave(q_tune_line, filename = here(plot_dir, "q_tune_line.png"), height = 5)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
