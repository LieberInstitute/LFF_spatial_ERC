## Louise Huuki-Myers, May 2026
## Explore Oligo + OPC in Xenium data

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")

library("tidyverse")
library("crumblr")
library("variancePartition")

data_dir <- here("processed-data", "21_Xenium", "17_xenium_OligoOPC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "17_xenium_OligoOPC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
