## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/f3538df2e932f537f8670bf708f2ff9434ef5d91/code/05_harmony_BayesSpace/05-BayesSpace_k_search.R

## Required libraries
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")
library("tidyverse")
library("HDF5Array")
library("getopt")
library("scater")
library("harmony")
library("here")
library("sessioninfo")


## Create output directories
dir_plots <- here("plots", "05_spe_correct_cluster", "05_BayesSpace_Marker")
if(!dir.exists(dir_plots)) dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_Marker")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
if(!dir.exists(here(data_dir,"clusters_BayesSpace"))) dir.create(here(data_dir, "clusters_BayesSpace"), showWarnings = FALSE, recursive = TRUE)


#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))

## Extract marker set

## pull spatialDLPFC marker genes, combine w/ top.svgs
dlpfc_layer_markers <- read_csv(here("processed-data", "00_project_prep", "06_marker_genes", "dlpfc_layers_top100.csv"))

marker_genes <- dlpfc_layer_markers |>
    filter(ensembl %in% rownames(spe), dataset == "spatialDLPFC") |>
    pull(ensembl) |>
    unique()

length(marker_genes) # [1] 752

#### Get reduced dims ####
message(Sys.time(), " - PCA")
spe <- scater::runPCA(spe, 
               subset_row = marker_genes,
               ncomponents = 50,
               name = paste0("PCA"))

#### Harmony Batch Correction ####

message(Sys.time(), " - Harmony")
spe <- harmony::RunHarmony(spe, "sample_id")

## export Marker PCA & Harmony 

marker_reduced_dims <- list(marker_pca = reducedDim(spe,"PCA"),
                            marker_harmony = reducedDim(spe,"HARMONY"))

# maybe use fwrite or other 
save(marker_reduced_dims, file = here(data_dir, "erc_Marker_reducedDims.Rdata") )


#### Cluster ####
## Set Seed
set.seed(20240905)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

## Fix colaname for row and col
spe$row <- spe$array_row
spe$col <- spe$array_col

walk(c(2, 9, 16), function(k){
    ## Run BayesSpace
    message(Sys.time(), " - Running spatialCluster: k=", k, ", dimred = ", dimred)
    spe <- BayesSpace::spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 20000)
    
    message(Sys.time(), " - Format and Export")
    
    table(spe$spatial.cluster)
    
    spe$bayesSpace_temp <- as.factor(spe$spatial.cluster)
    bayesSpace_name <- paste0("BayesSpace_",name,"_k", k_nice)
    colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name
    
    cluster_export(
        spe,
        bayesSpace_name,
        cluster_dir = here(data_dir, "clusters_BayesSpace")
    )
    
    ## Visualize BayesSpace results
    message(Sys.time(), " - Visualize clusters")
    sample_ids <- unique(spe$sample_id)
    
    ## create color pallet
    cols <- Polychrome::palette36.colors(k)
    if(k ==2) cols <- cols[seq(2)] ## fix return 3 colors bug
    
    names(cols) <- sort(unique(spe[[bayesSpace_name]]))
    
    #   Use 'vis_grid_clus' to preserve all spots (including overlaps)
    p_list = vis_grid_clus(
        spe = spe,
        clustervar = bayesSpace_name,
        pdf_file = here(dir_plots, paste0(bayesSpace_name, "-ALL.pdf")),
        sort_clust = FALSE,
        colors = cols,
        spatial = FALSE,
        point_size = 1.1
    )
    
    
})


# slurmjobs::job_single('05_BayesSpace_Marker', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace_Marker.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
