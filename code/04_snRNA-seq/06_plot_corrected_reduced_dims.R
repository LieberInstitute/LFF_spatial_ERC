library("SingleCellExperiment")
library("scater")
library("purrr")
library("HDF5Array")
library("viridis")
library("here")
library("sessioninfo")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

plot_dir <- here("plots", "04_snRNA-seq", "06_plot_corrected_reduced_dims")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

## Load HD5F
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "sce_harmony"))
sce

## swap rownames to gene names for plotting markers
rownames(sce) <- rowData(sce)$Symbol 


#### plot ####
## categorical
walk(c("sample_id", "seq_round", "exp_round","APOE"), ~my_plot_reduced_dim(sce, prefix = "sn", var_type = "cat", dimred = "UMAP", my_var = .x, sufix = "HARMONY"))
walk(c("sample_id", "seq_round", "exp_round","APOE"), ~my_plot_reduced_dim(sce, prefix = "sn", var_type = "cat", dimred = "TSNE", my_var = .x, sufix = "HARMONY"))

