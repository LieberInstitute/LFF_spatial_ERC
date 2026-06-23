## Louise Huuki-Myers, June 2026
## Run nearest neighbors to find Oligo.3-Astro microenvironments

#### Set Up ####

library("qs2")
library("SpatialExperiment")
library("tidyverse")
library("here")
library("sessioninfo")


data_dir <- here("processed-data", "21_Xenium", "20_xenium_Oligo3_Astro")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "20_xenium_Oligo3_Astro")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


#### Load data ####

message(Sys.time(), " - Load SPE data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
# filter to singlet cells
spe <- spe[,spe$spot_class == "singlet"]


