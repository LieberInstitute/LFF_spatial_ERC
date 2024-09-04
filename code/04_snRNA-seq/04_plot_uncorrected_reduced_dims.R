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
plot_dir <- here("plots", "04_snRNA-seq", "04_plot_uncorrected_reduced_dims")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

## Load "uncorrected" sce data (postQC, GLMPCA)
message(Sys.time(), " - load data")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_uncorrected"))

## swap rownames to gene names for plotting markers
rownames(sce) <- rowData(sce)$Symbol 

#### plot ####
## categorical
walk(c("sample_id", "seq_round", "exp_round","APOE"), ~my_plot_reduced_dim(sce, dimred = "UMAP", my_var = .x, sufix = "uncorrected"))
walk(c("sample_id", "seq_round", "exp_round","APOE"), ~my_plot_reduced_dim(sce, dimred = "TSNE", my_var = .x, sufix = "uncorrected"))

## continuous
walk(c("sum", "detected", "subsets_Mito_percent"), ~my_plot_reduced_dim(sce, dimred = "UMAP", cat_var = FALSE, my_var = .x, sufix = "uncorrected"))

walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(sce, dimred = "UMAP", cat_var = FALSE, my_var = .x, sufix = "expres_uncorrected"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
