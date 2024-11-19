## Louise Huuki-Myers, Nov 2024
## Run BayesSpace w/ spatialDLPFC marker genes - run clustering
## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/f3538df2e932f537f8670bf708f2ff9434ef5d91/code/05_harmony_BayesSpace/05-BayesSpace_k_search.R

## Required libraries
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")

## Choose k
k <- as.numeric(
    #   Only one of these environment variables will be defined, so grab the
    #   defined one (handle SGE or SLURM)
    paste0(Sys.getenv("SLURM_ARRAY_TASK_ID"), Sys.getenv("SGE_TASK_ID"))
)

# k <-9

k_nice <- sprintf("k%02d", k)
message("Run BayesSpace: ", k_nice)

## Create output directories
plot_dir <- here("plots", "05_spe_correct_cluster", "05_BayesSpace_Marker")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "05_BayesSpace_Marker")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
if(!dir.exists(here(data_dir,"clusters_BayesSpace"))) dir.create(here(data_dir, "clusters_BayesSpace"), showWarnings = FALSE, recursive = TRUE)


#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))

## load marker_reduced_dims
load(here(data_dir, "erc_Marker_reducedDims.Rdata"), verbose = TRUE)

reducedDim(spe, "HARMONY") <- marker_reduced_dims$marker_harmony
reducedDim(spe, "PCA") <- marker_reduced_dims$marker_pca

#### Cluster ####
message(Sys.time(), " - Set Up Clustering")
## Set Seed
set.seed(20240905)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

## Fix colaname for row and col
spe$row <- spe$array_row
spe$col <- spe$array_col

## Run BayesSpace
message(Sys.time(), " - Running spatialCluster: k=", k, ", dimred = HARMONY")
spe <- BayesSpace::spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 20000)

message(Sys.time(), " - Format and Export")

table(spe$spatial.cluster)

spe$bayesSpace_temp <- as.factor(spe$spatial.cluster)
bayesSpace_name <- paste0("BayesSpace_Markers_", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

if(!dir.exists(here(data_dir, bayesSpace_name))) dir.create(here(data_dir, bayesSpace_name))

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here(data_dir, bayesSpace_name)
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
    pdf_file = here(plot_dir, paste0(bayesSpace_name, "-ALL.pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 1.1
)

# slurmjobs::job_single('05_BayesSpace_Marker', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace_Marker.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
