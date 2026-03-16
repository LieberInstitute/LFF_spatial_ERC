## Louise Huuki-Myers, Jan 2026
## Correlate NMF factors with variables

#### Set up ####

library("SingleCellExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
# library("singlet")
library("ComplexHeatmap")
library("spatialLIBD")
library("RcppML")

data_dir <- here("processed-data", "20_NMF", "02_NMF_var_correlation")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "02_NMF_var_correlation")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### load data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## Drop Br1289
sce <- sce[, sce$BrNum != "Br1289"]

cell_type_colors <- metadata(sce)$cell_type_colors

# colors 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)


## NMF data
nmf <- readRDS(here("processed-data", "20_NMF", "01_runNMF", "ERC_sn_nmf.RDS"))

class(nmf)
# [1] "nmf"
# attr(,"package")
# [1] "RcppML"

#### Explore nmf ####
## 50 factors
dim(nmf$h)
# 50 116198

#### NMF summary ####

walk2(c("cell_type_broad", "cell_type_anno", "APOE", "APOE_carrier"),
      list(cell_type_colors$broad, cell_type_colors$anno,  APOE_genotype_colors, APOE_carrier_colors),
      function(c, color_values){
          cell_type_stats <- RcppML::summary(nmf, group_by = sce[[c]])
          
          summary_plot <- plot(cell_type_stats, stat = "sum") +
              scale_fill_manual(values = color_values)
          
          ggsave(summary_plot, filename = here(plot_dir, sprintf("nmf_summary-%s.png", c)), width = 10)
          
      })

## umap

set.seed(123)
umap <- data.frame(uwot::umap(t(nmf$h)))

head(umap)

umap$cell_type_anno <- sce$cell_type_anno
umap$APOE <- sce$APOE

nmf_UMAP <- ggplot(umap, aes(X1, X2, color = cell_type_anno)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$anno)

ggsave(nmf_UMAP, filename = here(plot_dir, "nmf_UMAP-cell_type_anno.png"), height = 10 , width = 12)

nmf_UMAP_apoe <- ggplot(umap, aes(X1, X2, color = APOE)) +
    geom_point(size = 1) +
    theme_bw() +
    scale_color_manual(values = APOE_genotype_colors)

ggsave(nmf_UMAP_apoe, filename = here(plot_dir, "nmf_UMAP-APOE.png"), height = 10 , width = 12)

#### Technical Vars ####
# data<-as.data.frame(sce$BrNum)
# colnames(data)<-'brain'
# onehot_brain <- reshape2::dcast(data = data, rownames(data) ~ brain, length)
# 
# head(onehot_brain)
# dim(onehot_brain)

tech_var_mod <- model.matrix(~Rin + exp_round + seq_round + sum + detected + subsets_Mito_percent, colData(sce))
dim(tech_var_mod)
head(tech_var_mod)

tech_var_mod <- as.data.frame(tech_var_mod)
tech_var_mod[,"(Intercept)"] <- NULL

colnames(tech_var_mod) <- gsub("roundround", "round", colnames(tech_var_mod))

nmf_tech_cor <- cor(t(nmf@h), tech_var_mod)

head(nmf_tech_cor)

## heatmap of technical variables
pdf(here(plot_dir, "nmf_tech_var_cor_heatmap.pdf"), height = 12)
Heatmap(nmf_tech_cor,
        name = "cor")
dev.off()

nmf_anno_tech_var <- annotate_registered_clusters(
    cor_stats_layer = abs(nmf_tech_cor),
    confidence_threshold = 0.4,
    cutoff_merge_ratio = 0.1
)

nmf_anno_tech_var |> filter(layer_confidence == "good")
#   cluster layer_confidence          layer_label
# 1    nmf2             good             detected
# 2    nmf7             good         sum/detected
# 3   nmf13             good         detected/sum
# 4   nmf15             good             detected
# 5   nmf17             good             detected
# 6   nmf25             good         sum/detected
# 7   nmf36             good           exp_round3
# 8   nmf45             good subsets_Mito_percent


#### Demographic Vars ####

demo_var_mod <- model.matrix(~APOE_carrier + Anc_Afr + Age + Sex, colData(sce))
head(demo_var_mod)

demo_var_mod <- as.data.frame(demo_var_mod)
demo_var_mod[,"(Intercept)"] <- NULL

nmf_demo_cor <- cor(t(nmf@h), demo_var_mod)

head(nmf_demo_cor)

## heatmap of demonical variables
pdf(here(plot_dir, "nmf_demo_var_cor_heatmap.pdf"), height = 12)
Heatmap(nmf_demo_cor)
dev.off()

nmf_anno_demo_var <- annotate_registered_clusters(
    cor_stats_layer = abs(nmf_demo_cor),
    confidence_threshold = 0.2,
    cutoff_merge_ratio = 0.1
)

nmf_anno_demo_var |> filter(layer_confidence == "good")
# cluster layer_confidence     layer_label
# 1    nmf5             good APOE_carrierE4+
# 2   nmf24             good            SexM
# 3   nmf26             good             Age
# 4   nmf36             good            SexM
# 5   nmf45             good            SexM
# 6   nmf46             good APOE_carrierE4+


#### BrNum ####
BrNum_var_mod <- model.matrix(~BrNum, colData(sce))
head(BrNum_var_mod)

BrNum_var_mod <- as.data.frame(BrNum_var_mod)
BrNum_var_mod[,"(Intercept)"] <- NULL

colnames(BrNum_var_mod) <- gsub("BrNumBr", "Br", colnames(BrNum_var_mod))

nmf_BrNum_cor <- cor(t(nmf@h), BrNum_var_mod)

head(nmf_BrNum_cor)

## demo annotations

sample_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

demo_anno <- sample_info |> 
    filter(BrNum %in% colnames(nmf_BrNum_cor)) |>
    select(BrNum, APOE, Ancestry, Sex, Age) |>
    column_to_rownames("BrNum")

demo_anno <- demo_anno[colnames(nmf_BrNum_cor), ]

demo_anno_col_ha <- HeatmapAnnotation(df = demo_anno,
                                      col = list(APOE = APOE_genotype_colors,
                                                 Sex = sex_colors,
                                                 Ancestry = ancestry_colors))

## heatmap of BrNum
pdf(here(plot_dir, "nmf_BrNum_var_cor_heatmap.pdf"), height = 12)
Heatmap(nmf_BrNum_cor,
        name = "cor",
        bottom_annotation = demo_anno_col_ha)
dev.off()

nmf_anno_BrNum_var <- annotate_registered_clusters(
    cor_stats_layer = nmf_BrNum_cor,
    confidence_threshold = 0.4,
    cutoff_merge_ratio = 0.1
)

nmf_anno_BrNum_var |> filter(layer_confidence == "good")
# cluster layer_confidence layer_label
# 1   nmf41             good      Br5599
# 2   nmf45             good      Br1691

#### Cell type ####

cellType_var_mod <- model.matrix(~cell_type_anno, colData(sce))
head(cellType_var_mod)

cellType_var_mod <- as.data.frame(cellType_var_mod)
cellType_var_mod[,"(Intercept)"] <- NULL

colnames(cellType_var_mod) <- gsub("cell_type_anno", "", colnames(cellType_var_mod))
levels(sce$cell_type_anno)[!levels(sce$cell_type_anno) %in% colnames(cellType_var_mod)]

cellType_var_mod$Astro.1 <- as.integer(sce$cell_type_anno == "Astro.1")

nmf_cellType_cor <- cor(t(nmf@h), cellType_var_mod)

head(nmf_cellType_cor)
dim(nmf_cellType_cor)

## cell type annotation
cell_type_colors <- metadata(sce)$cell_type_colors$anno

cellType_anno <- data.frame(cellType = colnames(nmf_cellType_cor))
cellType_col_ha <- HeatmapAnnotation(df = cellType_anno,
                                      col = list(cellType = cell_type_colors))

## heatmap of cellType variables
pdf(here(plot_dir, "nmf_cellType_var_cor_heatmap.pdf"), height = 12, width = 9)
Heatmap(nmf_cellType_cor,
        name = "cor",
        bottom_annotation = cellType_col_ha)
dev.off()

## best cell type for each nmf
nmf_anno_cellType_var <- annotate_registered_clusters(
    cor_stats_layer = nmf_cellType_cor,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
)

nmf_anno_cellType_var |> filter(layer_confidence == "good") |> arrange(layer_label)
nmf_anno_cellType_var |> filter(grepl("Oligo.3", layer_label))
# cluster layer_confidence layer_label
# 1   nmf46             poor    Oligo.3*

## flip and find top nmf for each cell type
nmf_anno_cellType_var2 <- annotate_registered_clusters(
    cor_stats_layer = t(nmf_cellType_cor),
    confidence_threshold = 0.4,
    cutoff_merge_ratio = 0.1
)

#### combine annotations ####

nmf_annotation <- data.frame(nmf = rownames(nmf@h))|> 
    left_join(nmf_anno_tech_var |> 
                  # filter(layer_confidence == "good") |> 
                  select(nmf = cluster, tech_var = layer_label)) |>
    left_join(nmf_anno_demo_var|> 
                  # filter(layer_confidence == "good") |> 
                  select(nmf = cluster, demo_var = layer_label))|>
    left_join(nmf_anno_BrNum_var|> 
                  # filter(layer_confidence == "good") |> 
                  select(nmf = cluster, BrNum = layer_label))|>
    left_join(nmf_anno_cellType_var |> 
                  # filter(layer_confidence == "good") |> 
                  select(nmf = cluster, cellType = layer_label))

write.csv(nmf_annotation, file = here(data_dir, "nmf_annotaions.csv"))

#### Find Marker ####

## make nmf sce object
sce_nmf <- SingleCellExperiment(colData = colData(sce), assays = list(logcounts = nmf$h))

nmf_ct_markers <- scran::findMarkers(sce_nmf, groups=sce_nmf$cell_type_anno)

nmf_ct_markers_tb <- map2_dfr(nmf_ct_markers, names(nmf_ct_markers), ~.x |>
             as.data.frame() |>
             rownames_to_column("Factor") |>
             transmute(Factor, Top, p.value, FDR, summary.logFC, cell_type = .y))

nmf_ct_markers_top <- nmf_ct_markers_tb |> group_by(cell_type) |> arrange(-summary.logFC) |> slice(1)

nmf_ct_markers_top |> print(n = 38)

nmf_anno_cellType_var |> filter(grepl("Astro.1", layer_label))



nmf_anno_cellType_var |>
    rename(Factor = cluster) |>
    left_join(nmf_ct_markers_top) |>
    arrange(-summary.logFC)

library(scDotPlot)

pdf(here(plot_dir, sprintf("nmf_cell_type_dotplot_ALL.pdf")), height = 12, width = 10)
scDotPlot(sce_nmf,
          features = rownames(sce_nmf),
          group = "cell_type_anno",
          scale = FALSE)

scDotPlot(sce_nmf,
          features = rownames(sce_nmf),
          group = "cell_type_anno",
          scale = TRUE)
dev.off()

nmf_anno_cellType_good <- nmf_anno_cellType_var2 |>
    filter(layer_confidence == "good"| cluster == "Oligo.3")|>
    separate(layer_label, into = c("top_factor", "other_factor"), sep = "/", extra = "merge") |>
    arrange(cluster) |>
    mutate(top_factor = gsub("\\*", "", top_factor))

select_ct_nmf <- unique(nmf_anno_cellType_good$top_factor)

pdf(here(plot_dir, sprintf("nmf_cell_type_dotplot_select.pdf")), height = 12, width = 10)
scDotPlot(sce_nmf,
          features = select_ct_nmf,
          group = "cell_type_anno",
          groupAnno = "cell_type_anno",
          annoColors = list("cell_type_anno" = cell_type_colors),
          scale = FALSE,
          clusterColumns = FALSE)

scDotPlot(sce_nmf,
          features = select_ct_nmf,
          group = "cell_type_anno",
          groupAnno = "cell_type_anno",
          annoColors = list("cell_type_anno" = cell_type_colors),
          scale = TRUE)
dev.off()



# slurmjobs::job_single('02_NMF_var_correlation', create_shell = TRUE, memory = '100G', command = "Rscript 02_NMF_var_correlation.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()






