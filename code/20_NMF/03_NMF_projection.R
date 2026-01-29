## Louise Huuki-Myers, Jan 2026
## Correlate NMF factors with variables

#### Set up ####

# library("SpatialExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("ComplexHeatmap")
library("spatialLIBD")

data_dir <- here("processed-data", "20_NMF", "03_NMF_projection")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "03_NMF_projection")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####

## snRNA-seq
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## spatail data
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

## Drop Br1289
sce <- sce[, sce$BrNum != "Br1289"]
spe <- spe[, sce$BrNum != "Br1289"]

#### NMF data ####
nmf <- readRDS(here("processed-data", "20_NMF", "01_runNMF", "ERC_sn_nmf.RDS"))

# extract patterns
patterns <- t(nmf@h)
colnames(patterns) <- paste("NMF", 1:75, sep = "_")

loadings <- x@w
rownames(loadings) <- rownames(sce)

