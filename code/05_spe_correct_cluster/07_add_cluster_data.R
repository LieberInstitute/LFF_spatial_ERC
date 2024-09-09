## September 2024, Louise Huuki-Myers
## Add clustering data to spe object 


## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("HDF5Array")
library("tidyverse")

## Load the data
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects","spe_postQC"))

#### Add missing colData ###
colData(spe)$APOE_carrier = ifelse(grepl("E2", spe$APOE), "E2+", "E4+")
table(spe$APOE_carrier)

## Load clustering data
cluster_fn <- list.files(here("processed-data", "05_spe_correct_cluster", "05_BayesSpace", "clusters_BayesSpace-PCA_Harmony", "clusters_BayesSpace"),
           recursive = TRUE,
           full.names = TRUE)


clusters <- map2(cluster_fn, map(str_split(cluster_fn, "/"), 12), function(file, name){
    sp = sprintf("Sp%02d",parse_number(name))
    clus <- read.csv(file) |>
        mutate(cluster = factor(paste0(sp, "D", str_pad(cluster, 2,  pad = "0"))))
    
    colnames(clus)[2] <- name
    return(clus)
})

clusters2 <- purrr::reduce(clusters, dplyr::left_join, by = c("key"))

## fix double BrNum in key
clusters2$key <- gsub("(^.*?_Br[0-9]+)_Br[0-9+]+$", "\\1", clusters2$key)
identical(colnames(spe), clusters2$key)
clusters2$key <- NULL

## bind tables
colData(spe) <- cbind(colData(spe), clusters2)

colnames(colData(spe))

#### Save data ####
message(Sys.time(), " - Saving HDF5 SPE")
saveHDF5SummarizedExperiment(
    spe,
    dir = here("processed-data", "spe_objects", "spe_ERC"),
    replace = TRUE
)

# slurmjobs::job_single('07_add_cluster_data', create_shell = TRUE, memory = '5G', command = "Rscript 07_add_cluster_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

