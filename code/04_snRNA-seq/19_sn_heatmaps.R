## Louise Huuki-Myers, March 2025
## Compare ERC spatial Domains to DLPFC layers

library("spatialLIBD")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "04_snRNA-seq", "19_sn_heatmaps")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "19_sn_heatmaps")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
## load psuedobulked sce data
sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk","sce_pseudobulk-cell_type_fine.rds"))
dim(sce_pb)

rownames(sce_pb) <- rowData(sce_pb)$Symbol

anno_notes <- readxl::read_excel(here("processed-data", "04_snRNA-seq", "13_sctype_final", "sctype_final-NOTES.xlsx"))

cell_anno_df <- colData(sce_pb) |>
    as.data.frame() |>
    select(sample_id, APOE, Ancestry, Sex, Age, cell_type_broad, cell_type_fine, cell_type_fine, ncells) |>
    left_join(anno_notes |> select(cell_type_fine, cell_type_anno = guess)) |>
    mutate(APOE = gsub("^(E[2,3,4])(E[2,3,4])","\\1/\\2", APOE))

cell_type_anno_levels <- cell_anno_df |> arrange(cell_type_broad, cell_type_anno) |> pull(cell_type_anno) |> unique()

rownames(cell_anno_df) <- colnames(sce_pb)

cell_anno_df$cell_type_anno <- factor(cell_anno_df$cell_type_anno, levels = cell_type_anno_levels)
cell_anno_df <- cell_anno_df |> arrange(cell_type_anno, sample_id)

sce_pb <- sce_pb[,rownames(cell_anno_df)]

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

setequal(names(APOE_genotype_colors), unique(cell_anno_df$APOE))

cell_type_colors_anno <- DeconvoBuddies::create_cell_colors(cell_types = cell_type_anno_levels,
                                                            pallet_name = "classic", 
                                                            split = "\\.")
cell_type_colors$anno <- cell_type_colors_anno$fine

## cell annotations simple
cell_row_ha_simple <- rowAnnotation(
    df = cell_anno_df |> select(cell_type_anno),
    col = list(cell_type_anno = cell_type_colors$anno)
)

cell_col_ha_simple <- HeatmapAnnotation(
    df = cell_anno_df |> select(cell_type_anno),
    col = list(cell_type_anno = cell_type_colors$anno)
)

# cell annotaions with APOE
cell_row_ha_APOE <- rowAnnotation(
    df = cell_anno_df |> select(cell_type_anno, APOE),
    col = list(cell_type_anno = cell_type_colors$anno,
               APOE = APOE_genotype_colors)
)

cell_col_ha_APOE <- HeatmapAnnotation(
    df = cell_anno_df |> select(cell_type_anno, APOE),
    col = list(cell_type_anno = cell_type_colors$anno,
               APOE = APOE_genotype_colors)
)

# cell annotations with all details
cell_row_ha_details <- rowAnnotation(
    df = cell_anno_df |> select(cell_type_anno, APOE, Ancestry, Sex, Age),
    col = list(cell_type_anno = cell_type_colors$anno,
               APOE = APOE_genotype_colors,
               Ancestry = ancestry_colors,
               Sex = sex_colors)
)

cell_col_ha_details <- HeatmapAnnotation(
    df = cell_anno_df |> select(cell_type_anno, APOE, Ancestry, Sex, Age),
    col = list(cell_type_anno = cell_type_colors$anno,
               APOE = APOE_genotype_colors,
               Ancestry = ancestry_colors,
               Sex = sex_colors)
)

#### Lit gene heatmap ####
## read in marker genes from lit
lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    mutate(in_data = gene_name %in% rowData(sce_pb)$Symbol,
           cell_type_broad = factor(cell_type_broad, levels = levels(sce_pb$cell_type_broad))) |>
    arrange(cell_type_broad)

## missing from pb data
lit_markers |> filter(is.na(cell_type_broad))
lit_markers |> filter(!in_data)
lit_markers |> filter(gene_name == "SYT1")

lit_markers <- lit_markers |> filter(in_data & !is.na(cell_type_broad))

nrow(lit_markers)

## gene annotation df

lit_gene_annotation <- lit_markers |>
    select(gene_name, cell_type, cell_type_broad) |>
    column_to_rownames("gene_name")

lit_gene_col_ha <- HeatmapAnnotation(df = lit_gene_annotation |> select(" " = cell_type_broad),
                                     col = list(" " = cell_type_colors$broad))

## extract z-scores
lit_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(lit_gene_annotation),]))
dim(lit_markers_zscore)
lit_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_sn_lit_gene_heatmap.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
                      name = "Z Score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      column_split = lit_gene_annotation$cell_type_broad,
                      # row_split = cell_anno_df$cell_type_anno,
                      right_annotation = cell_row_ha_simple,
                      bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_lit_gene_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        # column_split = lit_gene_annotation$cell_type_broad,
        # row_split = cell_anno_df$cell_type_anno,
        right_annotation = cell_row_ha_simple,
        bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_lit_gene_heatmap_cluster_details.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        # column_split = lit_gene_annotation$cell_type_broad,
        # row_split = cell_anno_df$cell_type_anno,
        right_annotation = cell_row_ha_details,
        bottom_annotation = lit_gene_col_ha
)
dev.off()


#### MeanRatio heatmap ####
## load marker gene data
load(here("processed-data", "04_snRNA-seq", "16_sn_MeanRatio", "MarkerStats_cell_type_fine.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(MeanRatio.rank <= 5, MeanRatio > 1, gene_name %in% rownames(sce_pb)) |>
    select(gene, MeanRatio.rank, cell_type_fine = cellType.target) |>
    left_join(anno_notes |> select(cell_type_fine, cell_type_anno = guess)) |>
    mutate(cell_type_anno = factor(cell_type_anno, levels = cell_type_anno_levels)) |>
    arrange(cell_type_anno) |>
    column_to_rownames("gene")

marker_stats_top |> count(cell_type_anno) |> arrange(n)


MR_gene_row_ha <- rowAnnotation(df = marker_stats_top |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors$anno))

MR_gene_col_ha <- HeatmapAnnotation(df = marker_stats_top |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors$anno))

## extract z-scores
MR_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(marker_stats_top),]))
dim(MR_markers_zscore)
MR_markers_zscore[1:5,1:5]

pdf(here(plot_dir, "ERC_sn_MR_gene_heatmap.pdf"), height = 14, width = 10)
Heatmap(t(MR_markers_zscore),
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$cell_type_anno,
        row_split = gsub("\\.","\n",marker_stats_top$cell_type_anno),
        bottom_annotation = cell_col_ha_details,
        right_annotation = MR_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_MR_gene_heatmap_cluster.pdf"), height = 12, width = 10)
Heatmap(t(MR_markers_zscore),
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$cell_type_anno,
        # row_split = cell_anno_df$cell_type_anno,
        bottom_annotation = cell_col_ha_simple,
        right_annotation = MR_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_MR_gene_heatmap_cluster_details.pdf"), height = 12, width = 10)
Heatmap(t(MR_markers_zscore),
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$cell_type_anno,
        # row_split = cell_anno_df$cell_type_anno,
        bottom_annotation = cell_col_ha_details,
        right_annotation = MR_gene_row_ha
)
dev.off()

#### Enrichment heatmap ####
## load marker gene data
modeling_results <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk", "modeling_results-cell_type_fine.rds"))

rowData(sce_pb)$gene_name <- rowData(sce_pb)$ID

enrichment_stats_top <- sig_genes_extract(
    n = 5,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pb
) |>
    select(ensembl, fdr, top, logFC, cell_type_fine = test) |>
    left_join(anno_notes |> select(cell_type_fine, cell_type_anno = guess)) |>
    mutate(cell_type_anno = factor(cell_type_anno, levels = cell_type_anno_levels)) |>
    arrange(cell_type_anno) |>
    column_to_rownames("ensembl")

enrichment_stats_top |> count(cell_type_anno) |> arrange(n)

## gene annotations
enrich_gene_row_ha <- rowAnnotation(df = enrichment_stats_top |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors$anno))

enrich_gene_col_ha <- HeatmapAnnotation(df = enrichment_stats_top |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors$anno))

## extract z-scores
enrich_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(enrichment_stats_top),]))
dim(enrich_markers_zscore)
enrich_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_sn_enrich_gene_heatmap.pdf"), height = 14, width = 10)
Heatmap(t(enrich_markers_zscore),
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$cell_type_anno,
        row_split = gsub("\\.","\n",enrichment_stats_top$cell_type_anno),
        bottom_annotation = cell_col_ha_details,
        right_annotation = enrich_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_enrich_gene_heatmap_cluster.pdf"), height = 12, width = 10)
Heatmap(t(enrich_markers_zscore),
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$cell_type_anno,
        # row_split = cell_anno_df$cell_type_anno,
        bottom_annotation = cell_col_ha_simple,
        right_annotation = enrich_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_enrich_gene_heatmap_cluster_details.pdf"), height = 12, width = 10)
Heatmap(t(enrich_markers_zscore),
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$cell_type_anno,
        # row_split = cell_anno_df$cell_type_anno,
        bottom_annotation = cell_col_ha_details,
        right_annotation = enrich_gene_row_ha
)
dev.off()


# slurmjobs::job_single('19_sn_heatmaps', create_shell = TRUE, memory = '25G', command = "Rcript 19_sn_heatmaps.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
