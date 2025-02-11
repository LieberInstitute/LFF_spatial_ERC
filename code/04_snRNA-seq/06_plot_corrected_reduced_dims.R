# Feb, 2025 - Louise Huuki-Myers
# Plot Harmony corrected reduced dimensions for snRNA-seq data

library("SingleCellExperiment")
library("scater")
library("purrr")
library("HDF5Array")
library("viridis")
library("here")
library("sessioninfo")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## set up plot dir
plot_dir <- here("plots", "04_snRNA-seq", "06_plot_corrected_reduced_dims")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

## Load "HARMONY" sce data (postQC, GLMPCA)
message(Sys.time(), " - load data")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_harmony"))

reducedDimNames(sce)

## swap rownames to gene names for plotting markers
rownames(sce) <- rowData(sce)$Symbol 

#### plot ####
## categorical
walk(c("sample_id", "seq_round", "exp_round","APOE"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", my_var = .x, var_type = "cat", sufix = "HARMONY"))
walk(c("sample_id", "seq_round", "exp_round","APOE"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "TSNE", my_var = .x, var_type = "cat", sufix = "HARMONY"))

## continuous
walk(c("sum", "detected", "subsets_Mito_percent"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", my_var = .x, var_type = "con", sufix = "HARMONY"))

walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", var_type = "express", my_var = .x, sufix = "HARMONY"))
walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(sce, prefix = "ERC_sn", dimred = "UMAP", var_type = "express", my_var = .x, sufix = "HARMONY"))

# slurmjobs::job_single('06_plot_corrected_reduced_dims', create_shell = TRUE, memory = '25G', command = "Rscript 06_plot_corrected_reduced_dims.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
