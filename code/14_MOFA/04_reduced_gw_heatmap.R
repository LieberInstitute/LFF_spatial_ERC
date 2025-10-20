library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

source(here("code", "13_compile_DGE", "logFC_heatmap.R"))

num_factors = 7
specific_factor = 'Factor3'
fdr_cutoff = 0.05
z_cutoff_genes = 1.96
max_genes = 2
my_views = c("Oligo.3", "Oligo.4", "Oligo.5", "Excit.L5.2", "WM_Sp09D06")

cell_colors_path = here(
    "processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"
)
spd_colors_path = here("processed-data", "SpD_colors.Rdata")
project_colors_path = here("processed-data", "project_colors.Rdata")
model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
sce_path = here(
    "processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn",
    "sce_pseudo_DGE-cell_type_anno.RDS"
)
spd_dge_path = here(
    "processed-data", "13_compile_DGE", "01_compile_DGE", "Visium",
    "DGE_results_carrier_Visium.Rds"
)
sn_dge_path = here(
    "processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine",
    "DGE_results_carrier_sn_fine.Rds"
)
plot_dir = here("plots", "14_MOFA", "04_reduced_gw_heatmap")

dir.create(file.path(plot_dir, 'GO'), recursive = TRUE, showWarnings = FALSE)

#   Prepare view colors
load(cell_colors_path, verbose = TRUE)
load(spd_colors_path, verbose = TRUE)
load(project_colors_path, verbose = TRUE)
names(SpD_colors) = sub('~', '_', names(SpD_colors))
view_colors = c(cell_type_colors$anno, SpD_colors)[my_views]
view_colors = c(view_colors, Multi = "grey30")

model = load_model(model_path)
sce = readRDS(sce_path)

gene_weights = get_geneweights(model = model, factor = specific_factor) |>
    left_join(
        rowData(sce) |>
            as.data.frame() |>
            dplyr::select(gene_id, gene_name),
        by = c("feature" = "gene_id")
    ) |>
    as_tibble() |>
    #   Z-score across genes within each view
    group_by(ctype) |>
    mutate(value = (value - mean(value)) / sd(value)) |>
    ungroup()

stopifnot(!any(is.na(gene_weights$gene_name)))

top_gene_weights = gene_weights |>
    filter(ctype %in% my_views, value > z_cutoff_genes)

stopifnot(nrow(top_gene_weights) > 0)

#   Only display up to max_genes per view
select_genes = top_gene_weights |>
    group_by(ctype) |>
    arrange(desc(value)) |>
    slice_head(n = max_genes) |>
    pull(feature)

################################################################################
#   Heatmap of top-weighted genes across all views
################################################################################

#   Annotation of rows-- each gene is labeled by the view it is a signature gene
#   for
view_table = top_gene_weights |>
    filter(feature %in% select_genes) |>
    group_by(gene_name) |> 
    summarise(n = n(), view = paste0(ctype, collapse = ", ")) |>
    mutate(
        view = factor(
            ifelse(n > 1, "Multi", view), levels = c("Multi", my_views)
        )
    ) |>
    dplyr::select(gene_name, view) |>
    column_to_rownames("gene_name") |>
    arrange(view)

view_table_row = rowAnnotation(df = view_table, col = list(view = view_colors))

## prep heatmap 
top_gw_value_matrix = gene_weights |>
    filter(feature %in% select_genes) |>
    dplyr::select(gene_name, ctype, value) |>
    pivot_wider(names_from = ctype, values_from = value) |>
    column_to_rownames("gene_name") |>
    as.matrix()
top_gw_value_matrix = top_gw_value_matrix[rownames(view_table), my_views]

min_val = min(top_gw_value_matrix, na.rm = TRUE)
max_val = max(top_gw_value_matrix, na.rm = TRUE)
max_abs_val = max(abs(min_val), abs(max_val))
col_fun = colorRamp2(
    c(-1 * max_abs_val, 0, max_abs_val),
    c("blue", "white", "red")
)

## weights for selected clusters only
pdf(
    file.path(plot_dir, sprintf('gene_weights_select_%s.pdf', specific_factor)),
    width = (length(my_views) / 4) + 2, height = nrow(view_table) / 4 + 1
)
Heatmap(
    top_gw_value_matrix,
    name = "feature\nweights\n(Z-score)",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = view_table_row,
    col = col_fun
)
dev.off()

#   Also produce a corresponding heatmap of logFC from DGE results
dge_data = bind_rows(readRDS(spd_dge_path), readRDS(sn_dge_path)) |>
    mutate(cluster = sub('~', '_', cluster)) |>
    filter(cluster %in% my_views)

cluster_levels = c(names(cell_type_colors$anno), names(SpD_colors))
logFC_Heatmap(
    dge_data, 
    gene_list = rownames(view_table),
    title = "Factor3", 
    h = nrow(view_table) / 4 + 1, 
    w = length(my_views) / 4 + 2, 
    cluster_col = FALSE,
    datatype = "MOFA",
    flip = TRUE,
    save = TRUE,
    order_genes = FALSE,
    row_anno = view_table_row
)

################################################################################
#   GO of top-weighted genes
################################################################################

do_go = function(gene_list, universe, plot_path) {
    for (ont_type in c("BP", "MF", "CC")) {
        go_obj = compareCluster(
            gene_list, fun = "enrichGO", universe = universe,
            OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
            pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE,
            keyType = "ENSEMBL"
        )

        if(!is.null(go_obj)) {
            go_obj@compareClusterResult = go_obj@compareClusterResult |>
                filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type)
            
            if(nrow(go_obj@compareClusterResult) > 0) {
                pdf(
                    sprintf(plot_path, ont_type),
                    height = as.integer(
                        min(11, nrow(go_obj@compareClusterResult) * 0.7)
                    ),
                )
                print(dotplot(go_obj, showCategory = 15))
                dev.off()
            }
        }
    }
}

#-------------------------------------------------------------------------------
#   One view at a time
#-------------------------------------------------------------------------------

for (this_view in my_views) {
    gene_list = list()
    gene_list[['up']] = top_gene_weights |>
        filter(value > 0, ctype == this_view) |>
        pull(feature)
    gene_list[['down']] = c()

    message(
        sprintf(
            'For view "%s", using %d up and %d down signature genes',
            this_view, length(gene_list[['up']]), length(gene_list[['down']])
        )
    )

    #   Only consider genes measured in this view
    universe = gene_weights |>
        filter(ctype == this_view) |>
        pull(feature)

    do_go(
        gene_list,
        universe,
        file.path(
            plot_dir, 'GO', sprintf('%%s_%s_%s.pdf', specific_factor, this_view)
        )
    )
}

#-------------------------------------------------------------------------------
#   Selected views pooled together
#-------------------------------------------------------------------------------

gene_list = list()
gene_list[['up']] = top_gene_weights |>
    filter(value > 0) |>
    pull(feature)
gene_list[['down']] = c()

#   Only consider genes measured in any of the selected views
universe = gene_weights |>
    filter(ctype %in% my_views) |>
    pull(feature) |>
    unique()

do_go(
    gene_list,
    universe,
    file.path(
        plot_dir, 'GO', sprintf('%%s_%s_together.pdf', specific_factor)
    )
)

session_info()
