## Louise Huuki-Myers, July 2025
## Process and plot RCTD spot deconvolution data

#### set up ####
library("ggplot2")
# library("SpatialExperiment")
# library("SummarizedExperiment")
library("spacexr")
library("HDF5Array")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "15_spot_deconvolution", "02_RCTD_results")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "15_spot_deconvolution", "02_RCTD_results")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load results ####
results_spe <- readRDS(here("processed-data", "15_spot_deconvolution", "01_runRCTD", "RCDT_est_prop.rds"))

# Assays have cell types as rows and pixels as columns


#### plot composition for each spot ####
plotAllWeights(results_spe, r = 0.05, lwd = 0, title = "Doublet Results")
