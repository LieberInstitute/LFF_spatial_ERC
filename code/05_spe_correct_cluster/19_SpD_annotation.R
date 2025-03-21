## December 2024, Louise Huuki-Myers
## Check cluster prevalence vs. groups in study

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("processed-data", "05_spe_correct_cluster", "19_SpD_annotation")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

table(spe$BayesSpace_SVGm_k09, spe$sample_id)

spe_pb_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k09.rds"))


as.data.frame(colData(spe)) |>
    group_by(BayesSpace_SVGm_k09) |>
    dplyr::summarise(n = n()) 

# DFplyr
# Error in `n()`:
#     ! Must only be used inside data-masking verbs like `mutate()`, `filter()`, and `group_by()`

cluster_prop <- as.data.frame(colData(spe)) |>
    group_by(sample_id, APOE, Ancestry, Sex, BayesSpace_SVGm_k09) |>
    summarise(n = n())  |>
    group_by(sample_id, APOE, Ancestry, Sex) |>
    mutate(cluster_prop = n/sum(n))


