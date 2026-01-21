## Louise Huuki-Myers, Jan 2026
## Correlate NMF factors with variables

#### Set up ####

library("SingleCellExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("singlet")
library("ComplexHeatmap")
library("spatialLIBD")

data_dir <- here("processed-data", "20_NMF", "02_NMF_var_correlation")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "02_NMF_var_correlation")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## Drop Br1289
sce <- sce[, sce$BrNum != "Br1289"]

## NMF data
nmf <- readRDS(here("processed-data", "20_NMF", "01_runNMF", "ERC_sn_nmf.RDS"))

class(nmf)

#### Explore nmf ####
## 66 factors
dim(nmf@h)
# 66 122004

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

colnames(tech_var_mod) <- gsub("BrNumBr", "Br", colnames(tech_var_mod))
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

#    cluster layer_confidence          layer_label
# 1     nmf1             good             detected
# 2     nmf4             good         sum/detected
# 3    nmf10             good         detected/sum
# 4    nmf12             good             detected
# 5    nmf17             good             detected
# 6    nmf19             good             detected
# 7    nmf30             good         sum/detected
# 8    nmf31             good         detected/sum
# 9    nmf37             good                 SexM
# 10   nmf51             good                 SexM
# 11   nmf53             good subsets_Mito_percent

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
    confidence_threshold = 0.4,
    cutoff_merge_ratio = 0.1
)

nmf_anno_demo_var |> filter(layer_confidence == "good")

#### BrNum ####

BrNum_var_mod <- model.matrix(~BrNum, colData(sce))
head(BrNum_var_mod)

BrNum_var_mod <- as.data.frame(BrNum_var_mod)
BrNum_var_mod[,"(Intercept)"] <- NULL

colnames(BrNum_var_mod) <- gsub("BrNumBr", "Br", colnames(BrNum_var_mod))

nmf_BrNum_cor <- cor(t(nmf@h), BrNum_var_mod)

head(nmf_BrNum_cor)

## demo annotations
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

sample_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

demo_anno <- sample_info |> 
    filter(BrNum %in% colnames(nmf_BrNum_cor)) |>
    select(BrNum, APOE, Ancestry, Sex, Age) |>
    column_to_rownames("BrNum")

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

#### Cell type ####

cellType_var_mod <- model.matrix(~cell_type_anno, colData(sce))
head(cellType_var_mod)

cellType_var_mod <- as.data.frame(cellType_var_mod)
cellType_var_mod[,"(Intercept)"] <- NULL

colnames(cellType_var_mod) <- gsub("cell_type_anno", "", colnames(cellType_var_mod))

nmf_cellType_cor <- cor(t(nmf@h), cellType_var_mod)

head(nmf_cellType_cor)
dim(nmf_cellType_cor)

## cell type annotation
cell_type_colors <- metadata(sce)$cell_type_colors$anno

cellType_anno <- data.frame(cellType = colnames(nmf_cellType_cor))
cellType_col_ha <- HeatmapAnnotation(df = cellType_anno,
                                      col = list(cellType = cell_type_colors))

## heatmap of cellType variables
pdf(here(plot_dir, "nmf_cellType_var_cor_heatmap.pdf"), height = 12, width = 8)
Heatmap(nmf_cellType_cor,
        name = "cor",
        bottom_annotation = cellType_col_ha)
dev.off()

nmf_anno_cellType_var <- annotate_registered_clusters(
    cor_stats_layer = nmf_cellType_cor,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
)

nmf_anno_cellType_var |> filter(layer_confidence == "good") |> arrange(layer_label)
nmf_anno_cellType_var |> filter(grepl("Oligo.3", layer_label))
# cluster layer_confidence     layer_label
# 1   nmf21             good Oligo.3/Oligo.4


#### combine annotations ####

nmf_annotation <- data.frame(nmf = rownames(nmf@h))|> 
    left_join(nmf_anno_tech_var |> 
                  filter(layer_confidence == "good") |> 
                  select(nmf = cluster, tech_var = layer_label)) |>
    left_join(nmf_anno_demo_var|> 
                  filter(layer_confidence == "good") |> 
                  select(nmf = cluster, demo_var = layer_label))|>
    left_join(nmf_anno_cellType_var |> 
                  filter(layer_confidence == "good") |> 
                  select(nmf = cluster, cellType = layer_label))




