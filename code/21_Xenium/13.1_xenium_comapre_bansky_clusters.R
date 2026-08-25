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
library("spatialLIBD")
library("scDotPlot")

#### define dirs ####
data_dir <- here("processed-data", "21_Xenium", "13.1_xenium_compare_bansky_clusters")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "13.1_xenium_compare_bansky_clusters")
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
rownames(spe) <- rowData(spe)$gene_id
## make APOE syntatic
spe$APOE_syn <- gsub("/", ".", spe$APOE)

modeling_results <- map(cluster_vars, ~registration_wrapper(
    sce = spe,
    var_registration = .x,
    var_sample_id = "sample_id",
    covars = c("APOE_syn", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    min_ncells = 10,
))

enrich_top <- map(modeling_results, ~sig_genes_extract(
    modeling_results = .x,
    sce_layer = spe,
    n = 10,
    model_type = "enrichment"
))

enrich_top_all <- map2_dfr(enrich_top, names(enrich_top), ~.x |> mutate(bansky = .y, .before = 1))

write_csv(enrich_top_all, file = here(data_dir, "banksy_clustering_enrich_top_k9-12.csv"))

enrich_top_wide <- map(enrich_top, ~.x |>
    group_by(cluster = test) |>
    summarise(top_enrich = paste(gene, collapse = ", ")))


vSpD_mod <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))

cor_layer <- map(modeling_results, ~layer_stat_cor(stats = .x$enrichment,
                                                   modeling_results = vSpD_mod,
                                                   model_type = "enrichment",
                                                   top_n = NULL)
)

anno <- map(cor_layer, ~annotate_registered_clusters(
    cor_stats_layer = .x,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.05 ## very strict annotation merge
))

anno_enrich <- map2(anno, enrich_top_wide, ~.x |> left_join(.y) )

map2(anno_enrich, names(anno_enrich), ~write_csv(.x, file = here(data_dir, sprintf("banksy_clustering_annotation_%s.csv", .y))))

pdf(here(plot_dir, "layer_stat_cor_x_vs_vSpD.pdf"))

walk(k_labels, ~print(layer_stat_cor_plot(cor_stats_layer = cor_layer[[.x]],
                                           cluster_rows =TRUE,
                                           cluster_columns = FALSE,
                                           annotation = anno[[.x]],
                                           query_colors = cluster_colors_by_k[[.x]]))
)
dev.off()


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

#### Layer expression plots ####

lit_markers <- read_csv(here("processed-data","05_spe_correct_cluster", "00_lit_marker_genes_layer", "lit_layer_marker_summary.csv")) |>
    filter(gene_name %in% rowData(spe)$gene_name) |>
    arrange(Layer)

rowData(spe)$LitMarker <- NULL
rowData(spe)$LitMarker <- lit_markers$Layer[match(rownames(spe), lit_markers$gene_name)] 
table(rowData(spe)$LitMarker)

layer_colors <- c(Vasc = '#FF99C0',
                  Layer1 = "#F0027F",
                  Layer2 = "#377EB8",
                  Layer3 = "#4DAF4A",
                  Layer4 = "#984EA3",
                  Layer5 = "#FFD700",
                  Layer5a = "#E9FF70",
                  Layer5b = "#FF9770",
                  Layer6 = "#FF7F00",
                  WM = "grey90",
                  GM = "grey25",
                  `Para-subiculum` = "brown",
                  Interneuron = "red")


pdf(here(plot_dir, "Xenium_test_SpD_dotplot_LitMarkers.pdf"))
map2(cluster_vars, cluster_colors_by_k, ~spe |>
    scDotPlot(features = lit_markers$gene_name,
              group = .x,
              groupAnno = .x,
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list(.x = .y,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE))
dev.off()

#### compare cell types ####

#### Add RCTD data ####
# 
message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

rctd_data@results$results_df <- rctd_data@results$results_df[colnames(spe),]

## enforce alignment - if any spe cell isn't in rctd_data, this indexing would
## produce NA-filled rows with mismatched row names, so this must hold before cbind
stopifnot(identical(colnames(spe), rownames(rctd_data@results$results_df)))

spe$cell_id <- colnames(spe)
colData(spe) <- cbind(colData(spe), rctd_data@results$results_df[colnames(spe),])

spe$cell_type_anno <- spe$first_type
spe$cell_type_broad <- factor(gsub("\\..*?$", "", spe$cell_type_anno), levels = c("Astro", "Macro", "Micro","Oligo", "OPC","Vasc","Excit","Inhib"))
table(spe$cell_type_broad)


library("ComplexHeatmap")

load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)

cell_v_SpX <- map(cluster_vars, ~table(spe[, spe$spot_class == "singlet"][[.x]], spe[, spe$spot_class == "singlet"]$cell_type_anno))

## proportion cell type (SpX rows sum to 1)
cell_prop_v_SpX <- map(cell_v_SpX, ~sweep(.x, 1, rowSums(.x), FUN = "/"))
rowSums(cell_prop_v_SpX[[1]])

## proportion SpX (cell type cols sum to 1)
cell_v_SpX_prop <- map(cell_prop_v_SpX, ~sweep(.x, 2, colSums(.x), FUN = "/"))
colSums(cell_v_SpX_prop[[1]])

## PLOT HEATMAPS
pdf(here(plot_dir, "Xenium_bansky_test_SpX_v_cell_type_heatmap.pdf"), width = 10)

## proportion cell type (SpX row sum to 1)
map(cell_prop_v_SpX, ~ComplexHeatmap::Heatmap(.x,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = TRUE,
                        cluster_rows = TRUE  ## cluster SpX
                        )
)

## proportion cell type (SpX row sum to 1)
map(cell_v_SpX_prop, ~ComplexHeatmap::Heatmap(.x,
                        name = "prop SpX\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = TRUE,
                        cluster_rows = TRUE  ## cluster SpX
                        )
)

dev.off()



# slurmjobs::job_single('13.1_xenium_compare_bansky_clusters', create_shell = TRUE, memory = '50G', command = "Rscript 13.1_xenium_compare_bansky_clusters.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
