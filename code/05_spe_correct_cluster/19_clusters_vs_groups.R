## December 2024, Louise Huuki-Myers
## Check cluster prevalence vs. groups in study

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("processed-data", "05_spe_correct_cluster", "19_clusters_vs_groups")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

spe_pb_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k09.rds"))

