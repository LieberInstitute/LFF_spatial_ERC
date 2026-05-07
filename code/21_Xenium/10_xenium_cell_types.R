## Louise Huuki-Myers, May 2026
## Explore RCTD results

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")

library("scran")
library("scater")

# library("tidyverse")
# library("crumblr")
# library("variancePartition")

data_dir <- here("processed-data", "21_Xenium", "11_xenium_crumblr")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "11_xenium_crumblr")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

names(rctd_data@results)

head(rctd_data@results$results_df)


#### Reduced dims ####

# use logcounts of your preferred normalization
# spe <- logNormCounts(spe, assay.type = "counts")

assay(spe, "logcounts") <- log2(assay(spe, "nucleus_normcounts") + 1)

message(Sys.time(), " - running PCA")
spe <- runPCA(spe,
              ncomponents = 30
)

message(Sys.time(), " - running TSNE")
spe <- runTSNE(spe, dimred = "PCA")

message(Sys.time(), " - running UMAP")
spe <- runUMAP(spe, dimred = "PCA")

# get variance explained
pca_var <- attr(reducedDim(spe, "PCA"), "varExplained")
pca_var_pct <- pca_var / sum(pca_var) * 100

# elbow plot
var_explained <- data.frame(PC = seq_along(pca_var_pct),
           var_explained = pca_var_pct) 

var_explained_elbow <- var_explained |>
    ggplot(aes(x = PC, y = var_explained)) +
    geom_point() +
    geom_line() +
    labs(x = "PC", y = "Variance Explained (%)") +
    theme_bw()


source(here("code", "utils", "my_plot_reduced_dim.R"))

#### plot ####
## categorical
walk(c("sample_id", "seq_round", "exp_round","APOE","quick_cluster"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", my_var = .x, var_type = "cat", sufix = "uncorrected"))
walk(c("sample_id", "seq_round", "exp_round","APOE","quick_cluster"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "TSNE", my_var = .x, var_type = "cat", sufix = "uncorrected"))

## continuous
walk(c("sum", "detected", "subsets_Mito_percent"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", my_var = .x, var_type = "con", sufix = "uncorrected"))

walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", var_type = "express", my_var = .x, sufix = "uncorrected"))
walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", var_type = "express", my_var = .x, sufix = "uncorrected"))


