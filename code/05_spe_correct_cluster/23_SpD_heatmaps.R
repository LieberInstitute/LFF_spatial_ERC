## Louise Huuki-Myers, March 2025
## Create Heatmaps of litature and data-driven marker genes for the spatial domains

library("spatialLIBD")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Set up dirs ####
# data_dir <- here("processed-data", "05_spe_correct_cluster", "23_SpD_heatmaps")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "23_SpD_heatmaps")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
## load psuedobulked spe data
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "25_SpD_pseudobulk","spe_pseudobulk_only-SpD_syn.rds"))
dim(spe_pb)

rownames(spe_pb) <- rowData(spe_pb)$gene_name

spe_pb_sample <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium","spe_pseudo_DGE.RDS"))
dim(spe_pb_sample)

rownames(spe_pb_sample) <- rowData(spe_pb_sample)$gene_name

spd_anno_df <- colData(spe_pb_sample) |>
    as.data.frame() |>
    select(sample_id, APOE, Ancestry, Sex, Age, SpD, ncells)

spd_levels <- levels(spe_pb_sample$SpD)

identical(rownames(spd_anno_df), colnames(spe_pb_sample))

spd_anno_df <- spd_anno_df |> arrange(SpD, sample_id)

spe_pb_sample <- spe_pb_sample[,rownames(spd_anno_df)]

## load colors
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

setequal(names(APOE_genotype_colors), unique(spd_anno_df$APOE))

#### create annotations ####
spd_row_ha_simple <- rowAnnotation(
    SpD = colData(spe_pb)$SpD,
    col = list(SpD = SpD_colors)
)

spd_col_ha_simple <- HeatmapAnnotation(
    SpD = colData(spe_pb)$SpD,
    col = list(SpD = SpD_colors)
)

# SpD + sample annotaions with APOE
spd_row_ha_APOE <- rowAnnotation(
    df = spd_anno_df |> select(SpD, APOE),
    col = list(SpD = SpD_colors,
               APOE = APOE_genotype_colors)
)

spd_col_ha_APOE <- HeatmapAnnotation(
    df = spd_anno_df |> select(SpD, APOE),
    col = list(SpD = SpD_colors,
               APOE = APOE_genotype_colors)
)

# cell annotations with all details
spd_row_ha_details <- rowAnnotation(
    df = spd_anno_df |> select(SpD, APOE, Ancestry, Sex, Age),
    col = list(SpD = SpD_colors,
               APOE = APOE_genotype_colors,
               Ancestry = ancestry_colors,
               Sex = sex_colors)
)

spd_col_ha_details <- HeatmapAnnotation(
    df = spd_anno_df |> select(SpD, APOE, Ancestry, Sex, Age),
    col = list(SpD = SpD_colors,
               APOE = APOE_genotype_colors,
               Ancestry = ancestry_colors,
               Sex = sex_colors)
)

#### Lit gene heatmap ####
message("** Lit genes **")
## read in marker genes from lit
lit_markers <- read_csv(here("processed-data","05_spe_correct_cluster", "00_lit_marker_genes_layer", "lit_layer_marker_summary.csv")) |>
    mutate(in_data = gene_name %in% rowData(spe_pb_sample)$gene_name) |>
    arrange(Layer)

## missing from pb data
lit_markers |> filter(!in_data)

lit_markers <- lit_markers |> filter(in_data)

nrow(lit_markers)

## gene annotation df

lit_gene_annotation <- lit_markers |>
    select(gene_name, Layer) |>
    column_to_rownames("gene_name")

unique(lit_gene_annotation$Layer)

## libd_layer_colors + more 
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
                  `Para-subiculum` = "brown")

lit_gene_col_ha <- HeatmapAnnotation(df = lit_gene_annotation |> select(" " = Layer),
                                     col = list(" " = layer_colors))


## extract z-scores
lit_markers_zscore <- scale(t(logcounts(spe_pb)[rownames(lit_gene_annotation),]))
dim(lit_markers_zscore)
lit_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_SpD_lit_gene_heatmap.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
                      name = "Z Score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = TRUE,
                      column_split = lit_gene_annotation$Layer,
                      # row_split = spd_anno_df$SpD,
                      right_annotation = spd_row_ha_simple,
                      bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_lit_gene_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        # column_split = lit_gene_annotation$Layer,
        # row_split = spd_anno_df$SpD,
        right_annotation = spd_row_ha_simple,
        bottom_annotation = lit_gene_col_ha
)
dev.off()

#### Lit gene heatmap - SpD + sample ####
## extract z-scores
lit_markers_zscore <- scale(t(logcounts(spe_pb_sample)[rownames(lit_gene_annotation),]))
dim(lit_markers_zscore)
lit_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_SpD_lit_gene_heatmap_sample.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
                      name = "Z Score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      # column_split = lit_gene_annotation$Layer,
                      # row_split = spd_anno_df$SpD,
                      right_annotation = spd_row_ha_APOE,
                      bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_lit_gene_heatmap_cluster_sample.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        # column_split = lit_gene_annotation$Layer,
        # row_split = spd_anno_df$SpD,
        right_annotation = spd_row_ha_APOE,
        bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_lit_gene_heatmap_cluster_details.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        # column_split = lit_gene_annotation$SpD,
        # row_split = spd_anno_df$SpD,
        right_annotation = spd_row_ha_details,
        bottom_annotation = lit_gene_col_ha
)
dev.off()

#### AD risk gene heatmap ####
## read in risk genes from OpenTargets data
AD_risk <- read_csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    filter(symbol %in% rowData(spe_pb)$gene_name) 

nrow(AD_risk)

## gene annotation df
AD_risk_annotation <- AD_risk |>
    select(symbol, eva) |>
    column_to_rownames("symbol") |>
    filter(eva > 0.5)

eva_colors = colorRamp2(c(0, 1), c("white", "darkcyan"))
AD_risk_col_ha <- HeatmapAnnotation(df = AD_risk_annotation, col = list(eva = eva_colors))

## extract z-scores - cell type pb
AD_risk_zscore <- scale(t(logcounts(spe_pb)[rownames(AD_risk_annotation),]))
dim(AD_risk_zscore)
AD_risk_zscore[1:5,1:5]

## plot heatmaps - cell type only
pdf(here(plot_dir, "ERC_SpD_AD_risk_heatmap.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        right_annotation = spd_row_ha_simple,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_AD_risk_heatmap_split.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = 2,
        row_split = 2,
        right_annotation = spd_row_ha_simple,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_AD_risk_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        right_annotation = spd_row_ha_simple,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

#### MeanRatio heatmap ####
message("** MeanRatio genes **")
## load marker gene data
load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_SpD.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(MeanRatio.rank <= 5, gene_name %in% rownames(spe_pb)) |>
    select(gene, MeanRatio.rank, SpD = cellType.target) |>
    mutate(SpD = factor(SpD, levels = spd_levels)) |>
    arrange(SpD) |>
    column_to_rownames("gene")

marker_stats_top |> count(SpD)


MR_gene_row_ha <- rowAnnotation(df = marker_stats_top |> select(" " = SpD),
                                     col = list(" " = SpD_colors))

MR_gene_col_ha <- HeatmapAnnotation(df = marker_stats_top |> select(" " = SpD),
                                     col = list(" " = SpD_colors))

## extract z-scores
MR_markers_zscore <- scale(t(logcounts(spe_pb)[rownames(marker_stats_top),]))
dim(MR_markers_zscore)
MR_markers_zscore[1:5,1:5]

pdf(here(plot_dir, "ERC_SpD_MR_gene_heatmap_sample.pdf"), height = 8, width = 10)
Heatmap(MR_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        # row_split = gsub("\\.","\n",marker_stats_top$SpD),
        right_annotation = spd_row_ha_simple,
        bottom_annotation = MR_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_MR_gene_heatmap_cluster_sample.pdf"), height = 8, width = 10)
Heatmap(MR_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        # row_split = gsub("\\.","\n",marker_stats_top$SpD),
        right_annotation = spd_row_ha_simple,
        bottom_annotation = MR_gene_col_ha
)
dev.off()

#### MeanRatio heatmap - SpD + sample ####
## extract z-scores
MR_markers_zscore <- scale(t(logcounts(spe_pb_sample)[rownames(marker_stats_top),]))
dim(MR_markers_zscore)
MR_markers_zscore[1:5,1:5]

pdf(here(plot_dir, "ERC_SpD_MR_gene_heatmap_sample.pdf"), height = 8, width = 10)
Heatmap(MR_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        # row_split = gsub("\\.","\n",marker_stats_top$SpD),
        right_annotation = spd_row_ha_APOE,
        bottom_annotation = MR_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_MR_gene_heatmap_cluster.pdf"), height = 8, width = 10)
Heatmap(MR_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        # row_split = gsub("\\.","\n",marker_stats_top$SpD),
        right_annotation = spd_row_ha_APOE,
        bottom_annotation = MR_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_MR_gene_heatmap_cluster_details.pdf"), height = 8, width = 10)
Heatmap(MR_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$SpD,
        # row_split = spd_anno_df$SpD,
        right_annotation = spd_row_ha_details,
        bottom_annotation = MR_gene_col_ha
)
dev.off()

#### Enrichment heatmap ####
message("** Enrichment genes **")
## load marker gene data
modeling_results <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))

rowData(spe_pb_sample)$gene_name <- rowData(spe_pb_sample)$ID

enrichment_stats_top <- sig_genes_extract(
    n = 5,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    spe_layer = spe_pb
) |>
    select(ensembl, fdr, top, logFC, SpD = test) |>
    mutate(SpD = factor(gsub("_", "~", SpD), levels = spd_levels)) |>
    arrange(SpD) |>
    column_to_rownames("ensembl")

enrichment_stats_top |> count(SpD) |> arrange(n)

## gene annotations
enrich_gene_row_ha <- rowAnnotation(df = enrichment_stats_top |> select(" " = SpD),
                                     col = list(" " = SpD_colors))

enrich_gene_col_ha <- HeatmapAnnotation(df = enrichment_stats_top |> select(" " = SpD),
                                     col = list(" " = SpD_colors))

## extract z-scores
enrich_markers_zscore <- scale(t(logcounts(spe_pb)[rownames(enrichment_stats_top),]))
dim(enrich_markers_zscore)
enrich_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_SpD_enrich_gene_heatmap.pdf"), height = 8, width = 10)
Heatmap(enrich_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        right_annotation = spd_row_ha_simple,
        bottom_annotation = enrich_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_enrich_gene_heatmap_cluster.pdf"), height = 8, width = 10)
Heatmap(enrich_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        right_annotation = spd_row_ha_simple,
        bottom_annotation = enrich_gene_col_ha
)
dev.off()

#### Enrichment Heatmap - SpD + Sample ####
## extract z-scores
enrich_markers_zscore <- scale(t(logcounts(spe_pb_sample)[rownames(enrichment_stats_top),]))
dim(enrich_markers_zscore)
enrich_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_SpD_enrich_gene_heatmap_sample.pdf"), height = 8, width = 10)
Heatmap(enrich_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        right_annotation = spd_row_ha_APOE,
        bottom_annotation = enrich_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_enrich_gene_heatmap_cluster_sample.pdf"), height = 8, width = 10)
Heatmap(enrich_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$SpD,
        right_annotation = spd_row_ha_APOE,
        bottom_annotation = enrich_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_SpD_enrich_gene_heatmap_cluster_details.pdf"), height = 8, width = 10)
Heatmap(enrich_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$SpD,
        # row_split = spd_anno_df$SpD,
        right_annotation = spd_row_ha_details,
        bottom_annotation = enrich_gene_col_ha
)
dev.off()


# slurmjobs::job_single('23_SpD_heatmaps', create_shell = TRUE, memory = '25G', command = "Rscript 23_SpD_heatmaps.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
