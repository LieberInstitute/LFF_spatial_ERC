## Louise Huuki-Myers, Aug 2025
## MOFA gene weight heatmap

#### set up ####
library("SpatialExperiment")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test 
# opt$datatype = "sn_broad"

data_dir <- here("processed-data", "14_MOFA", "02_MOFA_gw_heatmaps", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "14_MOFA", "02_MOFA_gw_heatmaps", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## colors
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
