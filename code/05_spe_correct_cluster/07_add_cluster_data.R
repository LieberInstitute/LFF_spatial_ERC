## September 2024, Louise Huuki-Myers (Update Nov w/ SVGm data)
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
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_PCA_SVG"))

#### Add missing colData ###
colData(spe)$APOE_carrier = ifelse(grepl("E2", spe$APOE), "E2+", "E4+")
table(spe$APOE_carrier)

#### Load BayesSpace clustering data ####
bayes_cluster_fn <- list.files(here("processed-data", "05_spe_correct_cluster", "05_BayesSpace", "clusters_BayesSpace-SVGm", "clusters_BayesSpace"),
           recursive = TRUE,
           full.names = TRUE)

bayes_clusters <- map2(bayes_cluster_fn, map(str_split(bayes_cluster_fn, "/"), 12), function(file, name){
    sp = sprintf("Sp%02d",parse_number(name))
    clus <- read.csv(file) |>
        mutate(cluster = factor(paste0(sp, "D", str_pad(cluster, 2,  pad = "0"))))
    
    colnames(clus)[2] <- name
    return(clus)
})

bayes_cluster_markers_fn <- list.files(here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_Marker"),
                                       recursive = TRUE,
                                       full.names = TRUE,
                                       pattern = "clusters.csv")
## only select k02 and k11
bayes_cluster_markers_fn <- bayes_cluster_markers_fn[grepl("k02|k11", bayes_cluster_markers_fn)]

bayes_clusters_markers <- map2(bayes_cluster_markers_fn, map(str_split(bayes_cluster_markers_fn, "/"), 11), function(file, name){
    sp = sprintf("Sp%02d",parse_number(name))
    clus <- read.csv(file) |>
        mutate(cluster = factor(paste0(sp, "D", str_pad(cluster, 2,  pad = "0"))))
    
    colnames(clus)[2] <- name
    return(clus)
})


bayes_clusters_tab <- purrr::reduce(c(bayes_clusters, bayes_clusters_markers), dplyr::left_join, by = c("key"))
head(bayes_clusters_tab)

## fix double BrNum in key
bayes_clusters_tab$key <- gsub("(^.*?_Br[0-9]+)_Br[0-9+]+$", "\\1", bayes_clusters_tab$key)
identical(colnames(spe), bayes_clusters_tab$key)
bayes_clusters_tab$key <- NULL

bayes_clusters_tab[1:5,1:5]


#### PRECAST DATA ####
# precast_tab <-read.csv(file = here("processed-data", "05_spe_correct_cluster", "06_PRECAST", "PRECAST_bayes_clusters.csv"), row.names = 1)
# identical(rownames(precast_tab), colnames(spe))
# 
# precast_tab[1:5, 1:5]

## only add BayesSpace data for now
all_clusters <- cbind(bayes_clusters_tab)
# all_clusters <- cbind(bayes_clusters_tab, precast_tab)

dim(all_clusters)

## bind tables
colData(spe) <- cbind(colData(spe), all_clusters)

colnames(colData(spe))


#### Save data ####
message(Sys.time(), " - Saving HDF5 SPE")
saveHDF5SummarizedExperiment(
    spe,
    dir = here("processed-data", "spe_objects", "spe_ERC"),
    replace = TRUE
)

# slurmjobs::job_single('07_add_cluster_data', create_shell = TRUE, memory = '10G', command = "Rscript 07_add_cluster_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

