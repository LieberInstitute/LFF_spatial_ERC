## October 2024, Louise Huuki-Myers
## plot layer markers

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "05_spe_correct_cluster", "16_plot_layer_markers")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "16_plot_layer_markers")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

rownames(spe) <- rowData(spe)$gene_name
