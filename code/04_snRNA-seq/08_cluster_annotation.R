## Louise Huuki-Myers, September 2024
## Add cluster data to sce object, explore and annotate clusters

library("SingleCellExperiment")
library("jaffelab")
# library("scater")
# library("scran")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## Prep directories
plot_dir <- here("processed-data", "04_snRNA-seq", "08_cluster_annotation")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "08_cluster_annotation")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "sce_harmony"))
sce

## load cluster outputs
cluster_fn <- list.files(here("processed-data", "04_snRNA-seq", "07_cluster_sn"), full.names = TRUE)
names(cluster_fn) <- ss(basename(cluster_fn), "_", 3)

clusters <- map(cluster_fn, ~get(load(.x)))

map(clusters, max)
map(clusters, table)

sce$snn_k20 <- sprintf("k20c%02d", clusters$k20)
table(sce$snn_k20)

walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn", var_type = "cat", dimred = .x, my_var = "snn_k20", sufix = "HARMONY"))


#### iSEE export ####
lobstr::obj_size(sce)
# 292.63 MB

assays(sce)$binomial_deviance_residuals <- NULL

saveRDS(sce, here("code", "06_iSEE_app", "sce_ERC_iSEE.rds"))



