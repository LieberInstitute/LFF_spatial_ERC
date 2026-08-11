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
bayes_cluster_fn <- list.files(here("processed-data", "05_spe_correct_cluster", "05_BayesSpace", "clusters_BayesSpace-SVGm_200k", "clusters_BayesSpace"),
           recursive = TRUE,
           full.names = TRUE)

map(bayes_cluster_fn, ~as.Date(file.info(.x)$mtime))

bayes_clusters <- map2(bayes_cluster_fn, map(str_split(bayes_cluster_fn, "/"), 12), function(file, name){
    name <- gsub("_200k", "", name)
    sp = sprintf("Sp%02d",parse_number(name))
    clus <- read.csv(file) |>
        mutate(cluster = factor(paste0(sp, "D", str_pad(cluster, 2,  pad = "0"))))
    
    colnames(clus)[2] <- name
    return(clus)
})


bayes_clusters_tab <- purrr::reduce(bayes_clusters, dplyr::left_join, by = c("key"))
head(bayes_clusters_tab)

## fix double BrNum in key
bayes_clusters_tab$key <- gsub("(^.*?_Br[0-9]+)_Br[0-9+]+$", "\\1", bayes_clusters_tab$key)
identical(colnames(spe), bayes_clusters_tab$key)
bayes_clusters_tab$key <- NULL

bayes_clusters_tab[1:5,1:5]
dim(bayes_clusters_tab)
# [1] 122202     29

#### Add annotations ####
if(FALSE){
    ## TODO add to metadata
    load(here("processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC", "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"), verbose = TRUE)
    # cor_anno$BayesSpace_SVGm_k09$layer_anno$HumanPilot
    
    cluster_anno_k9 <- cor_anno$BayesSpace_SVGm_k09$layer_anno$HumanPilot |>
        mutate(BayesSpace_SVGm_k09_anno = paste0(layer_label_simple, "~", cluster))
    
    bayes_clusters_tab$BayesSpace_SVGm_k09_anno <- cluster_anno_k9$BayesSpace_SVGm_k09_anno[match(bayes_clusters_tab$BayesSpace_SVGm_k09, cluster_anno_k9$cluster)]
}

## bind tables
colData(spe) <- cbind(colData(spe), bayes_clusters_tab)

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

