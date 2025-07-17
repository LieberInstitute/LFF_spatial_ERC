## Louise Huuki-Myers, June 2025
## Create vis gene + vis clus plots of all samples

#### Set up ####

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

## if we need QC annotations
# message(Sys.time(), " - Load QC HDF5 SPE")
# spe_qc <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
# spe_qc

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

#### Supp fig plots ####

sample_info <- colData(spe) |>
    as.data.frame() |>
    select(sample_id, APOE) |>
    unique() |>
    mutate(APOE_syn = gsub("/", "", APOE)) |>
    as_tibble()

samples_by_apoe <- map(rafalib::splitit(sample_info$APOE_syn), ~sample_info$sample_id[.x])

walk(names(samples_by_apoe), function(apoe){
    
    samples <- samples_by_apoe[[apoe]]
    
    ## set low margin
    no_margin <- margin(t = 0,  # Top margin
                        r = 0,  # Right margin
                        b = 0,  # Bottom margin
                        l = 0) # Left margin
    
    ## H&E Column
    # test <- spatialLIBD::vis_image(spe, sampleid = "Br1556")
    
    image_col <- Reduce("/", map(samples, ~spatialLIBD::vis_image(spe, sampleid = .x) + 
                                     coord_fixed() +
                                     labs(y = .x, title = NULL) +
                                     theme(plot.margin = no_margin)
    )
    )
    
    # ggsave(image_col, filename = here(plot_dir, "test_image_col.png"), width = 3, height = 11)
    
    ## Sum UMI Column
    UMI_col <- Reduce("/", map(samples, ~vis_gene(spe, sampleid = .x, 
                                                  point_size = 1, 
                                                  geneid = "sum_umi") + 
                                   theme(plot.margin = no_margin) +
                                   labs(title = NULL)
    )
    )
    
    # ggsave(UMI_col, filename = here(plot_dir, "test_UMI_col.png"), width = 6, height = 18)
    
    
    ## RBG Column
    RGB_col <- Reduce("/", map(samples, ~vis_RGB(spe, sampleid = .x, 
                                                 point_size = 1) + 
                                   theme(legend.position = "None",
                                         plot.margin = no_margin) +
                                   labs(title = NULL)
    )
    )
    
    # ggsave(RGB_col, filename = here(plot_dir, "test_RGB_col.png"), width = 4, height = 18)
    
    
    ## SpD Column
    SpD_col <- Reduce("/", map(samples, ~vis_clus(spe, sampleid = .x, 
                                                  point_size = 1, 
                                                  clustervar = "SpD", 
                                                  colors = SpD_colors) + 
                                   theme(legend.position = "None",
                                         plot.margin = no_margin) +
                                   labs(title = NULL)
    )
    )
    
    # ggsave(SpD_col, filename = here(plot_dir, "test_SpD_col.png"), width = 3, height = 11)
    
    ## combine cols into grid
    sample_grid <- Reduce("|", list(image_col, UMI_col, RGB_col, SpD_col))
    
    nsamples <- length(samples)
    
    ggsave(sample_grid, 
           filename = here(plot_dir, sprintf("Visium_sample_grid_%s.png", apoe)), 
           width = nsamples*7/3, 
           height = nsamples*3)
    
    ggsave(sample_grid, 
           filename = here(plot_dir, sprintf("Visium_sample_grid_%s.pdf", apoe)), 
           width = nsamples*7/3, 
           height = nsamples*3)
    
    
})


