## Louise Huuki-Myers, Jan 2026
## Explore existing pannel, build custom list


#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")

data_dir <- here("processed-data", "21_Xenium", "01_xenium_panel")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "01_xenium_panel")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


