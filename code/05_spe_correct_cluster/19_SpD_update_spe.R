## Louise Huuki-Myers, March 2025
## Add SpD annotaions to spe, explore and annotate SpDs

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

plot_dir <- here("processed-data", "05_spe_correct_cluster", "19_SpD_update_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "19_SpD_update_spe")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))
spe

table(spe$BayesSpace_SVGm_k09, spe$sample_id)

## clean up colData
colnames(colData(spe))

colData(spe) <- colData(spe)[,!grepl("PRECAST", colnames(colData(sce)))]

#### load SpD ####

spd_anno <- readxl::read_excel(here("processed-data","05_spe_correct_cluster", "10_spatial_registration_DLPFC", "ERC_SpD_spatial_registration_Annotations.xlsx"))

spd_anno <- spd_anno |>
    mutate(SpD = fct_reorder(paste0(Annotation, "~", cluster), order))


#### Save data ####
message(Sys.time(), " - Saving HDF5 SPE")
saveHDF5SummarizedExperiment(
    spe,
    dir = here("processed-data", "spe_objects", "spe_ERC"),
    replace = TRUE
)

