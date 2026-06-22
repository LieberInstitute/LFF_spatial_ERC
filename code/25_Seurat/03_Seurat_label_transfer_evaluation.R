## Louise Huuki-Myers, Jun 2026
## Evaluate label transfer between the snRNA-seq data and the seurat objects for each donor
## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/145f74fb821473fb36d6845c5f6331ea0024f089/code/24_xenium/08_snRNA_xenium_transfer.R

#### Set Up ####

library("sessioninfo")
library("Seurat")
library("here")
library("qs2")
library("bluster")
library("ComplexHeatmap")
library("SpatialExperiment")
library("tidyverse")

data_dir <- here("processed-data", "25_Seurat", "03_Seurat_label_transfer_evaluation")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "25_Seurat", "03_Seurat_label_transfer_evaluation")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 

#### Load Data ####
message(Sys.time(), "- Load xenium data")
prediction <- readRDS(here("processed-data", "25_Seurat", "02_Seurat_label_transfer", "xenium_snRNA_predictions.Rds"))

table(prediction$predicted.id)

message(Sys.time(), " - Load QS2 xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
spe <- spe[,spe$spot_class == "singlet"]

identical(rownames(prediction), colnames(spe))

## Add Seurat data
spe$Seurat_cell_type <- factor(prediction$predicted.id, levels = levels(spe$cell_type_anno))
spe$Seurat_score_max <- prediction$prediction.score.max

#### Agreement with RCDT ####
table(spe$cell_type_anno == spe$Seurat_cell_type)
# FALSE   TRUE 
# 31937 270342

sum(spe$cell_type_anno == spe$Seurat_cell_type)/ncol(spe) # 89% agreement on cell type calls
sum(spe$Seurat_score_max >= 0.5)/ncol(spe) # 92% max score > 0.5

table(spe$cell_type_anno, spe$Seurat_cell_type)

jacc.mat <- linkClustersMatrix(spe$cell_type_anno, spe$Seurat_cell_type)

identical(rownames(jacc.mat), colnames(jacc.mat))

ct_broad <- gsub("\\..*","", rownames(jacc.mat))

## column barplot
col_barplot = HeatmapAnnotation(seurat_n = anno_barplot(as.numeric(table(spe$Seurat_cell_type)),
                                                        gp = gpar(fill = cell_type_colors$anno)))

row_barplot = rowAnnotation(rctd_n = anno_barplot(as.numeric(table(spe$cell_type_anno)),
                                                  gp = gpar(fill = cell_type_colors$anno)))


pdf(here(plot_dir, "xenium_jacc_matrix_RCTD_v_Seurat.pdf"), height = 10, width = 10)
Heatmap(jacc.mat,
        name = "Correspondence",
        col = c("black", viridisLite::plasma(100)),
        na_col = "black",
        row_split = ct_broad,
        column_split = ct_broad,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = col_barplot,
        right_annotation = row_barplot,
        row_title = "RCTD predictions",
        column_title = "Seurat predictions"
)
dev.off()

#### RCTD output ####
# message(Sys.time(), "- Load rctd data")
# rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))
# 
# names(rctd_data@results)
# 
# head(rctd_data@results$results_df)
# 
# summary(rctd_data@results$results_df$min_score)
# 
# summary(spe$singlet_score)

#### Score boxplots ####

prediction$RCTD_cell_type <- spe$cell_type_anno
prediction$RCTD_cell_type_broad <- spe$cell_type_broad

prediction <- prediction |>
    mutate(Seurat_cell_type = factor(prediction$predicted.id, levels = levels(spe$cell_type_anno)),
           Seurat_cell_type_broad = factor(gsub("\\..*", "", prediction$predicted.id), levels = levels(spe$cell_type_broad)),
           match_cell_type_anno = RCTD_cell_type == Seurat_cell_type,
           match_cell_type_broad = RCTD_cell_type_broad == Seurat_cell_type_broad
           )

max_score_box <- ggplot(prediction,
                        aes(x = Seurat_cell_type, 
                            y = prediction.score.max,
                            fill = Seurat_cell_type)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.50) +
    theme_bw() +
    facet_wrap(~Seurat_cell_type_broad, scales = "free_x") +
    scale_fill_manual(values = cell_type_colors$anno) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = "Predicted Cell Type",
         y = "Max Score")


ggsave(max_score_box,filename = here(plot_dir, "Seurat_transfer_Max_score_boxplot.png"))


max_score_box_match <- ggplot(prediction,
                        aes(x = RCTD_cell_type, 
                            y = prediction.score.max,
                            fill = RCTD_cell_type, 
                            color = match_cell_type_anno)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.50) +
    theme_bw() +
    facet_wrap(~RCTD_cell_type_broad, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cell_type_colors$anno) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") 

ggsave(max_score_box_match, filename = here(plot_dir, "Seurat_transfer_Max_score_boxplot_RCTD_match.png"), width = 12)


max_score_density <- ggplot(prediction,
                        aes(x = prediction.score.max,
                            fill = RCTD_cell_type, 
                            color = match_cell_type_anno)) +
    geom_density() +
    geom_hline(yintercept = 0.50) +
    theme_bw() +
    facet_wrap(~RCTD_cell_type) +
    scale_fill_manual(values = cell_type_colors$anno) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") 

ggsave(max_score_density, filename = here(plot_dir, "Seurat_transfer_Max_score_density_RCTD_match.png"), width = 12)

#### SEURAT pca ####
head(rownames(xen))
head(rownames(xen))

test <- FeaturePlot(xen, features = c("MOBP"))
