## Compare Banksy clusters (k09-k12) across resolutions, Xenium data
## Adapted from 18_compare_clusters.R (Visium SVGm resolution comparison),
## applied to the clust_HARMONY_kmeans{9,10,11,12} clusters produced by
## 13_xenium_bansky_embedding.R

library("tidyverse")
library("bluster")
library("ComplexHeatmap")
library("qs2")
library("here")
library("sessioninfo")
library("SpatialExperiment")

#### define dirs ####
data_dir <- here("processed-data", "21_Xenium", "13.1_xenium_comapre_bansky_clusters")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "13.1_xenium_comapre_bansky_clusters")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load prelim data (from 13_xenium_bansky_embedding.R) ####
message(Sys.time(), " - Load prelim spe data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "spe_xenium_bansky_prelim.qs2"))

k_values <- 9:12
k_labels <- sprintf("k%02d", k_values)

cluster_vars <- sprintf("clust_HARMONY_kmeans%d", k_values)
names(cluster_vars) <- k_labels

## sanity check the expected clustering columns are present
stopifnot(all(cluster_vars %in% colnames(colData(spe))))

#### Generic per-cluster colors ####
cluster_colors <- c("#f62062",
                    "#f45e28",
                    "#cf9800",
                    "#608d00",
                    "#01e090",
                    "#3ed9e6",
                    "#0064ca",
                    "#9215a3",
                    "#ff8ee2",
                    "black",
                    "grey",
                    "brown")

make_generic_colors <- function(k, palette = cluster_colors) {
    pad_width <- nchar(as.character(k))
    setNames(palette[seq_len(k)], sprintf("xSp%dD%0*d", k, pad_width, seq_len(k)))
}
cluster_colors_by_k <- map(k_values, make_generic_colors)
names(cluster_colors_by_k) <- k_labels

# $k09
# xSp9D1    xSp9D2    xSp9D3    xSp9D4    xSp9D5    xSp9D6    xSp9D7    xSp9D8    xSp9D9
# "#f62062" "#f45e28" "#cf9800" "#608d00" "#01e090" "#3ed9e6" "#0064ca" "#9215a3" "#ff8ee2"

#### Model pseudobulk ####


#### Compare resolutions with a Jaccard-style correspondence matrix ####
k_pairs <- combn(k_labels, 2, simplify = FALSE)
names(k_pairs) <- map_chr(k_pairs, ~paste(.x, collapse = "_"))

jacc_mat_list <- map(k_pairs, function(pair) {
    linkClustersMatrix(colData(spe)[[cluster_vars[[pair[1]]]]], colData(spe)[[cluster_vars[[pair[2]]]]])
})

jacc_mat_long <- map(jacc_mat_list, ~ .x |>
    reshape2::melt() |>
    dplyr::rename(cluster_row = Var1, cluster_col = Var2, Jacc = value))

## for each cluster in the first resolution, its single best-matching cluster
## in the second resolution
best_match <- imap_dfr(jacc_mat_long, ~ .x |>
    group_by(cluster_row) |>
    slice_max(Jacc) |>
    ungroup() )

write_csv(best_match, file = here(data_dir, "bansky_best_jaccard_match.csv"))

#### Jaccard heatmaps, annotated with each resolution's own cluster colors ####

#' Build a row (or column) annotation for a Jaccard heatmap, named after the
#' resolution it represents and colored from that resolution's
#' cluster_colors_by_k entry, so k09/k10/k11/k12 each keep consistent colors.
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
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        bottom_annotation = col_ha,
        col = c("black", viridisLite::plasma(100)),
        na_col = "black",
        row_title = row_k,
        column_title = col_k
    )
}

iwalk(jacc_mat_list, function(jacc_mat, pair_name) {
    row_k <- k_pairs[[pair_name]][1]
    col_k <- k_pairs[[pair_name]][2]

    pdf(here(plot_dir, sprintf("xenium_bansky_jacc_mat_%s.pdf", pair_name)))
    print(plot_jacc_heatmap(jacc_mat, row_k, col_k, cluster_colors_by_k))
    dev.off()
})

# slurmjobs::job_single('13.1_xenium_comapre_bansky_clusters', create_shell = TRUE, memory = '50G', command = "Rscript 13.1_xenium_comapre_bansky_clusters.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
