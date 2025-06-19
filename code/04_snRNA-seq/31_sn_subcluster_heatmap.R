## Louise Huuki-Myers, May 2025
## Compare ERC spatial Domains to DLPFC layers

library("spatialLIBD")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")
library("circlize")

#### Set up dirs ####
data_dir <- here("processed-data", "04_snRNA-seq", "31_sn_subcluster_heatmap")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "31_sn_subcluster_heatmap")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
## load psuedobulked (by cluster) sce data
sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_ERC_subcluster_pseudobulk_only-cell_type_anno.rds"))
dim(sce_pb)

rownames(sce_pb) <- rowData(sce_pb)$gene_name

# Broad cell types
sce_pb_broad <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_ERC_subcluster_pseudobulk_only-cell_type_broad.rds"))
rownames(sce_pb_broad) <- rowData(sce_pb_broad)$gene_name

## load psuedobulked (by cluster + sample) sce data
sce_pb_sample <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk","sce_subcluster_pseudobulk-cell_type_anno.rds"))
dim(sce_pb_sample)

rownames(sce_pb_sample) <- rowData(sce_pb_sample)$gene_name

table(sce_pb$cell_type_anno)

#### Prep cell type annotations ####
## Pull annotations 
cell_anno_df <- colData(sce_pb_sample) |>
    as.data.frame() |>
    select(sample_id, APOE, Ancestry, Sex, Age, cell_type_broad, cell_type_anno, ncells) |>
    mutate(APOE = gsub("^(E[2,3,4])(E[2,3,4])","\\1/\\2", APOE),
           cell_type_anno = droplevels(cell_type_anno)) |> 
    arrange(cell_type_anno, sample_id)

all(rownames(cell_anno_df) == colnames(sce_pb_sample))

cell_type_anno_levels <- levels(sce_pb$cell_type_anno)

# sce_pb <- sce_pb[,rownames(cell_anno_df)]

## load colors
cell_type_colors <- metadata(sce_pb)$cell_type_colors

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

## Check colors
cell_type_colors$anno <- cell_type_colors$anno[cell_type_anno_levels]
setequal(names(cell_type_colors$anno), cell_type_anno_levels)
# setdiff(names(cell_type_colors$anno), cell_type_anno_levels)

setequal(names(APOE_genotype_colors), unique(cell_anno_df$APOE))

## cell + sample annotations simple
cell_row_ha_simple <- rowAnnotation(
    cell_type_anno = colData(sce_pb)$cell_type_anno,
    col = list(cell_type_anno = cell_type_colors$anno)
)

cell_col_ha_simple <- HeatmapAnnotation(
    cell_type_anno = colData(sce_pb)$cell_type_anno,
    col = list(cell_type_anno = cell_type_colors$anno)
)
## cell + sample annotations broad
cell_row_ha_broad <- rowAnnotation(
    cell_type_broad = colData(sce_pb_broad)$cell_type_broad,
    col = list(cell_type_broad = cell_type_colors$broad)
)

cell_col_ha_broad <- HeatmapAnnotation(
    cell_type_broad = colData(sce_pb_broad)$cell_type_anno,
    col = list(cell_type_broad = cell_type_colors$anno)
)

## cell + sample annotations simple
cell_sample_row_ha_simple <- rowAnnotation(
    df = cell_anno_df |> select(cell_type_anno),
    col = list(cell_type_anno = cell_type_colors$anno)
)

cell_sample_col_ha_simple <- HeatmapAnnotation(
    df = cell_anno_df |> select(cell_type_anno),
    col = list(cell_type_anno = cell_type_colors$anno)
)

# cell + sample annotations with APOE
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
    mutate(in_data = gene_name %in% rowData(sce_pb)$gene_name,
           cell_type_broad = factor(cell_type_broad, levels = levels(sce_pb$cell_type_broad))) |>
    arrange(cell_type_broad)

## missing from pb data
lit_markers |> filter(is.na(cell_type_broad))
lit_markers |> filter(!in_data)

lit_markers <- lit_markers |> filter(in_data & !is.na(cell_type_broad))

lit_markers |> dplyr::count(cell_type_broad)

nrow(lit_markers)

## gene annotation df
lit_gene_annotation <- lit_markers |>
    select(gene_name, cell_type, cell_type_broad) |>
    column_to_rownames("gene_name")

lit_gene_col_ha <- HeatmapAnnotation(df = lit_gene_annotation |> select(" " = cell_type_broad),
                                     col = list(" " = cell_type_colors$broad))

## extract z-scores - cell type pb
lit_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(lit_gene_annotation),]))
dim(lit_markers_zscore)
lit_markers_zscore[1:5,1:5]

summary(lit_markers_zscore[,2])

## plot heatmaps - cell type only
pdf(here(plot_dir, "ERC_sn_subcluster_lit_gene_heatmap.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
                      name = "Z Score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      column_split = lit_gene_annotation$cell_type_broad,
                      # row_split = cell_anno_df$cell_type_anno,
                      right_annotation = cell_row_ha_simple,
                      bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_lit_gene_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        # column_split = lit_gene_annotation$cell_type_broad,
        # row_split = cell_anno_df$cell_type_anno,
        right_annotation = cell_row_ha_simple,
        bottom_annotation = lit_gene_col_ha
)
dev.off()


#### Lit heatmaps - cell type + sample ####
## subset gene annotatoins 
lit_gene_annotation <- lit_gene_annotation[rownames(lit_gene_annotation) %in% rownames(sce_pb_sample),]
nrow(lit_gene_annotation)

lit_gene_col_ha <- HeatmapAnnotation(df = lit_gene_annotation |> select(" " = cell_type_broad),
                                     col = list(" " = cell_type_colors$broad))

## extract z-scores - cell type + sample pb
lit_markers_zscore_sample <- scale(t(logcounts(sce_pb_sample)[rownames(lit_gene_annotation),]))
dim(lit_markers_zscore_sample)

## plot heatmaps
pdf(here(plot_dir, "ERC_sn_subcluster_lit_gene_heatmap_sample.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore_sample,
                      name = "Z Score",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      # column_split = lit_gene_annotation$cell_type_broad,
                      # row_split = cell_anno_df$cell_type_anno,
                      right_annotation = cell_sample_row_ha_simple,
                      bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_lit_gene_heatmap_cluster_sample.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore_sample,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        # column_split = lit_gene_annotation$cell_type_broad,
        # row_split = cell_anno_df$cell_type_anno,
        right_annotation = cell_sample_row_ha_simple,
        bottom_annotation = lit_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_lit_gene_heatmap_cluster_details.pdf"), height = 8, width = 11)
Heatmap(lit_markers_zscore_sample,
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

#### AD risk gene heatmap ####
## read in risk genes from OpenTargets data
AD_risk <- read_csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    filter(symbol %in% rowData(sce_pb)$gene_name) 

nrow(AD_risk)

## gene annotation df
AD_risk_annotation <- AD_risk |>
    select(symbol, eva) |>
    column_to_rownames("symbol")

eva_colors = colorRamp2(c(0, 1), c("white", "darkcyan"))
AD_risk_col_ha <- HeatmapAnnotation(df = AD_risk_annotation, col = list(eva = eva_colors))

## extract z-scores - cell type pb
AD_risk_zscore <- scale(t(logcounts(sce_pb)[rownames(AD_risk_annotation),]))
dim(AD_risk_zscore)
AD_risk_zscore[1:5,1:5]

## save risk gene sub clusters
write.csv(AD_risk_zscore, file = here(data_dir, "AD_risk_zscore_sn_subcluster.csv"), row.names = FALSE)

risk_gene_max_cell_type <- AD_risk_zscore |>
    as.data.frame() |>
    rownames_to_column("cell_type") |>
    pivot_longer(!cell_type, names_to = "gene", values_to = "zscore") |>
    group_by(gene) |>
    slice_max(zscore)

summary(risk_gene_max_cell_type$zscore)

spatialLIBD::annotate_registered_clusters(
    t(AD_risk_zscore),
    confidence_threshold = 2,
    cutoff_merge_ratio = 0.1
)

## plot heatmaps - cell type only
pdf(here(plot_dir, "ERC_sn_subcluster_AD_risk_heatmap.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        right_annotation = cell_row_ha_simple,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_AD_risk_heatmap_split.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = 5,
        row_split = 4,
        right_annotation = cell_row_ha_simple,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_AD_risk_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        right_annotation = cell_row_ha_simple,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

#### AD risk gene heatmap - broad ####

sce_pb_broad

## extract z-scores - cell type pb
AD_risk_zscore_broad <- scale(t(logcounts(sce_pb_broad)[rownames(AD_risk_annotation),]))
dim(AD_risk_zscore_broad)
AD_risk_zscore_broad[1:5,1:5]

## save risk gene sub clusters
write.csv(AD_risk_zscore_broad, file = here(data_dir, "AD_risk_zscore_sn_broad.csv"), row.names = FALSE)

risk_gene_max_cell_type <- AD_risk_zscore_broad |>
    as.data.frame() |>
    rownames_to_column("cell_type") |>
    pivot_longer(!cell_type, names_to = "gene", values_to = "zscore") |>
    group_by(gene) |>
    slice_max(zscore)

summary(risk_gene_max_cell_type$zscore)

risk_gene_cell_type_association <- spatialLIBD::annotate_registered_clusters(
    t(AD_risk_zscore_broad),
    confidence_threshold = 1,
    cutoff_merge_ratio = 0.1
)

write.csv(risk_gene_cell_type_association, file = here(data_dir, "risk_gene_cell_type_association.csv"))

## plot heatmaps - cell type broad
pdf(here(plot_dir, "ERC_sn_broad_AD_risk_heatmap.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore_broad,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        right_annotation = cell_row_ha_broad,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_broad_AD_risk_heatmap_split.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore_broad,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = 3,
        row_split = 3,
        right_annotation = cell_row_ha_broad,
        bottom_annotation = AD_risk_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_broad_AD_risk_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(AD_risk_zscore_broad,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        right_annotation = cell_row_ha_broad,
        bottom_annotation = AD_risk_col_ha
)
dev.off()


#### Excit subtype heatmap  ####
lit_markers_Excit_subtype <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes","Davila-Velderrain_Excit_subtype_markers.csv"))

subtype_gene_annotation <- lit_markers_Excit_subtype |>
    mutate(in_data = gene_name %in% rowData(sce_pb)$gene_name) |>
    filter(in_data, is.na(note) | note != "exclude") |>
    select(gene_name, cell_type) |>
    column_to_rownames("gene_name") |>
    arrange(cell_type)

subtype_gene_annotation |> count(cell_type)

subtype_gene_col_ha <- HeatmapAnnotation(df = subtype_gene_annotation |> select(" " = cell_type),
                                         col = list(" " = c(Excit_L2 = "#377EB8",
                                                                  Excit_L2_3 = "aquamarine",
                                                                  Excit_L3 = "#4DAF4A",
                                                                  Excit_L3b = "darkgreen",
                                                                  Excit_L5 = "#FFD700",
                                                                  # Excit_L5b = "#FF9770",
                                                                  Excit_L6 = "#FF7F00")))

## Excit only cell type annotations
cell_excit_row_ha_simple <- rowAnnotation(
    cell_type_anno = sce_pb$cell_type_anno[grepl("Excit",sce_pb$cell_type_anno)],
    col = list(cell_type_anno = cell_type_colors$anno[grepl("Excit", names(cell_type_colors$anno))])
)

cell_excit_row_ha_simple_sample <- rowAnnotation(
    df = cell_anno_df |> filter(grepl("Excit", cell_type_anno)) |> select(cell_type_anno),
    col = list(cell_type_anno = cell_type_colors$anno[grepl("Excit", names(cell_type_colors$anno))])
)

## extract z-scores
subtype_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(subtype_gene_annotation), grep("Excit", colnames(sce_pb))]))
dim(subtype_markers_zscore)
subtype_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_sn_subcluster_Excit_subtype_gene_heatmap.pdf"), height = 8, width = 11)
Heatmap(subtype_markers_zscore,
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        # row_split = cell_anno_df |> filter(grepl("Excit", cell_type_anno)) |> select(cell_type_anno),
        right_annotation = cell_excit_row_ha_simple,
        bottom_annotation = subtype_gene_col_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_Excit_subtype_gene_heatmap_cluster.pdf"), height = 8, width = 11)
Heatmap(subtype_markers_zscore,
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        # row_split = cell_anno_df |> filter(grepl("Excit", cell_type_anno)) |> select(cell_type_anno),
        right_annotation = cell_excit_row_ha_simple,
        bottom_annotation = subtype_gene_col_ha
)
dev.off()



#### MeanRatio heatmap ####
## load marker gene data
load(here("processed-data", "04_snRNA-seq", "16_sn_MeanRatio", "MarkerStats_cell_type_anno.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(MeanRatio.rank <= 5, MeanRatio > 1, gene_name %in% rownames(sce_pb)) |>
    select(gene_name, MeanRatio.rank, cell_type_anno = cellType.target) |>
    mutate(cell_type_anno = factor(cell_type_anno, levels = cell_type_anno_levels)) |>
    arrange(cell_type_anno) |>
    column_to_rownames("gene_name")

marker_stats_top |> count(cell_type_anno) |> arrange(n)


MR_gene_row_ha <- rowAnnotation(df = marker_stats_top |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors$anno))

MR_gene_col_ha <- HeatmapAnnotation(df = marker_stats_top |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors$anno))

## extract z-scores
MR_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(marker_stats_top),]))
dim(MR_markers_zscore)
MR_markers_zscore[1:5,1:5]

pdf(here(plot_dir, "ERC_sn_subcluster_MR_gene_heatmap.pdf"), height = 14, width = 10)
Heatmap(t(MR_markers_zscore),
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$cell_type_anno,
        row_split = gsub("\\.","\n",marker_stats_top$cell_type_anno),
        # bottom_annotation = MR_gene_col_ha, # flip
        bottom_annotation = cell_col_ha_simple,
        # right_annotation = cell_row_ha_simple #flip
        right_annotation = MR_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_MR_gene_heatmap_cluster.pdf"), height = 12, width = 10)
Heatmap(t(MR_markers_zscore),
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_split = gsub("\\.","\n",marker_stats_top$cell_type_anno),
        bottom_annotation = cell_col_ha_simple,
        right_annotation = MR_gene_row_ha
)
dev.off()

#### Mean Ratio heatmap - cell type + sample ####
## extract z-scores
MR_markers_zscore_sample <- scale(t(logcounts(sce_pb_sample)[rownames(marker_stats_top),]))
dim(MR_markers_zscore_sample)
MR_markers_zscore_sample[1:5,1:5]

pdf(here(plot_dir, "ERC_sn_subcluster_MR_gene_heatmap_sample.pdf"), height = 14, width = 10)
Heatmap(t(MR_markers_zscore_sample),
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

pdf(here(plot_dir, "ERC_sn_subcluster_MR_gene_heatmap_cluster_sample.pdf"), height = 12, width = 10)
Heatmap(t(MR_markers_zscore),
        name = "Z Score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 9),
        # column_split = marker_stats_top$cell_type_anno,
        # row_split = cell_anno_df$cell_type_anno,
        bottom_annotation = cell_sample_col_ha_simple,
        right_annotation = MR_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_MR_gene_heatmap_cluster_details.pdf"), height = 12, width = 10)
Heatmap(t(MR_markers_zscore_sample),
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
modeling_results <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_anno.rds"))

enrichment_stats_top <- sig_genes_extract(
    n = 5,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pb
) |>
    select(ensembl, fdr, top, logFC, cell_type_anno = test) |>
    arrange(cell_type_anno) 

enrichment_stats_top |> count(cell_type_anno) |> arrange(n)

## write out top genes
write.csv(enrichment_stats_top, file = here(data_dir, "sce_subcluster_enrichment_top5.csv"), row.names = FALSE)

## some genes have multiple cell types - simplify
enrichment_stats_top_unique <- enrichment_stats_top |>
    group_by(ensembl) |>
    mutate(cell_type_broad = jaffelab::ss(cell_type_anno, "\\.")) |>
    summarise(n = n(),
              ct = list(cell_type_anno),
              n_ctb = length(unique(cell_type_broad)), 
              ctb = paste(unique(cell_type_broad), collapse = "-")
              ) |>
    mutate(cell_type_anno = ifelse(n == 1, unlist(ct), ctb)) |>
    ungroup() |>
    arrange(cell_type_anno) |>
    column_to_rownames("ensembl")

enrichment_stats_top_unique |> count(cell_type_anno) |> arrange(n)

cell_type_colors_enrich <- c(cell_type_colors$broad, cell_type_colors$anno)[unique(enrichment_stats_top_unique$cell_type_anno)]

## gene annotations
enrich_gene_row_ha <- rowAnnotation(df = enrichment_stats_top_unique |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors_enrich))

enrich_gene_col_ha <- HeatmapAnnotation(df = enrichment_stats_top_unique |> select(" " = cell_type_anno),
                                     col = list(" " = cell_type_colors_enrich))

## extract z-scores
enrich_markers_zscore <- scale(t(logcounts(sce_pb)[rownames(enrichment_stats_top_unique),]))
dim(enrich_markers_zscore)
enrich_markers_zscore[1:5,1:5]

## plot heatmaps
pdf(here(plot_dir, "ERC_sn_subcluster_enrich_gene_heatmap.pdf"), height = 16, width = 10)
Heatmap(t(enrich_markers_zscore),
        name = "Z Score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = grid::gpar(fontsize = 8),
        # row_split_gp = grid::gpar(fontsize = 10),
        # column_split = marker_stats_top$cell_type_anno,
        # row_split = gsub("\\.","\n",enrichment_stats_top_unique$cell_type_anno),
        bottom_annotation = cell_col_ha_simple,
        right_annotation = enrich_gene_row_ha
)
dev.off()

pdf(here(plot_dir, "ERC_sn_subcluster_enrich_gene_heatmap_cluster.pdf"), height = 12, width = 10)
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


## TODO do by sample

pdf(here(plot_dir, "ERC_sn_subcluster_enrich_gene_heatmap_cluster_details.pdf"), height = 12, width = 10)
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


# slurmjobs::job_single('31_sn_subcluster_heatmap', create_shell = TRUE, memory = '10G', command = "Rscript 31_sn_subcluster_heatmap.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
