## Louise Huuki-Myers, March 2025
## Create Main fig plots from SPE object

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("cowplot")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

plot_dir <- here("plots", "05_spe_correct_cluster", "22_SpD_clean_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

table(spe$SpD, spe$sample_id)

#### load SpD annotations ####
spd_anno <- readxl::read_excel(here("processed-data","05_spe_correct_cluster", "10_spatial_registration_DLPFC", "ERC_SpD_spatial_registration_Annotations.xlsx"))

spd_anno$cluster %in% spe$BayesSpace_SVGm_k09

anno_table <- spd_anno |>
    mutate(Annotation = fct_reorder(Annotation, order)) |>
    mutate(SpD = fct_reorder(paste0(Annotation, "~", cluster), order)) |>
    select(cluster, Annotation, SpD) 

#### Define colors for SpD ####

SpD_colors <- c("Vasc~Sp09D08" = "#E05AD2",
                "L1~Sp09D05" = "#0220DE",
                "L2.3~Sp09D01" = "#FEAF16",
                "L3~Sp09D02" = "#00BCF9",
                "L4.inhib~Sp09D09" = "#C82100",
                "L5~Sp09D03" = "#16FF32",
                "L6~Sp09D04" = "#116A52",
                "WM.uf~Sp09D07" = "#E4E1E3",
                "WM~Sp09D06" = "#581009")


#### Plot reduced dims ####
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(spe,
                          prefix = "ERC_spe",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "SpD",
                          color_pal = SpD_colors))


# slurmjobs::job_single('22_SpD_clean_plots', create_shell = TRUE, memory = '25G', command = "Rscript 22_SpD_clean_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
