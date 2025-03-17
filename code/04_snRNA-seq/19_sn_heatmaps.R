## Louise Huuki-Myers, March 2025
## Compare ERC spatial Domains to DLPFC layers

library("spatialLIBD")
library("purrr")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "04_snRNA-seq", "18_sn_spatial_registration")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "18_sn_spatial_registration")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load psuedobulked sce data
sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk","sce_pseudobulk-cell_type_fine.rds"))
dim(sce_pb)

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)


# slurmjobs::job_single('18_sn_spatial_registration', create_shell = TRUE, memory = '25G', command = "Rcript 18_sn_spatial_registration.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
