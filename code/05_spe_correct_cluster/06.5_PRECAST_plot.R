## October 2024, Louise Huuki-Myers
## Plot spot plots for PRECAST output

library("HDF5Array")
library("spatialLIBD")
library("tidyverse")
        
#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "06_PRECAST")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "06.5_PRECAST_plot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

## PRECAST data
precast_output <- list.files(data_dir)