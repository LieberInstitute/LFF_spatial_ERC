## Louise Huuki-Myers, June 2025
## Create vis gene + vis clus plots of all samples

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("patchwork")

## source reduced dims function
source(here("code", "05_spe_correct_cluster", "vis_RGB.R"))

## set up dirs
plot_dir <- here("plots", "05_spe_correct_cluster", "29_SpD_vis_grid")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "05_spe_correct_cluster", "29_SpD_vis_grid")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

#### poster plots ####
poster_samples <- c("Br1556","Br5415", "Br1039")

# spatialLIBD::vis_image(spe, sampleid = "Br1556")

## H&E image
image_row <- Reduce("+", map(poster_samples, ~spatialLIBD::vis_image(spe, sampleid = .x) + coord_fixed()))
ggsave(image_row, filename = here(plot_dir, "poster_image_row.png"), width = 9, height = 4)

## marker gene RGB
rgb_row <-  Reduce("+", map(poster_samples, ~vis_RGB(spe, sampleid = .x, point_size = 1) + labs(title = NULL)))
ggsave(rgb_row, filename = here(plot_dir, "poster_rgb_row.png"), width = 9, height = 4)

rgb_legned <- vis_RGB_legend_simple()
ggsave(rgb_legned, filename = here(plot_dir, "rgb_legend_simple.png"))


## SpDs
SpD_row <- Reduce("+",
                  map(poster_samples, 
                      ~vis_clus(spe, sampleid = .x, 
                                point_size = 1, 
                                clustervar = "SpD", 
                                colors = SpD_colors) + 
                          theme(legend.position = "None") +
                          labs(title = NULL))
)

ggsave(SpD_row, filename = here(plot_dir, "poster_SpD_row.png"), width = 9, height = 4)

## grid
poster_grid <- image_row/rgb_row/SpD_row
ggsave(,filename = here(plot_dir, "poster_grid.png"), height = 10, width = 10)

