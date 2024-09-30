## September 2024, Louise Huuki-Myers
## run BayesSpace clustering with 10x prelim clusters 
## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/f3538df2e932f537f8670bf708f2ff9434ef5d91/code/05_harmony_BayesSpace/05-BayesSpace_k_search.R

## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")
library("tidyverse")
library("HDF5Array")
library("getopt")

## get input file from params
spec <- matrix(
    c(
        c("spe", "dimred", "name", "prelim"),
        c("s", "d", "n", "p"),
        c("1", "2", "3", "4"),
        c("character", "character", "character", "character"),
        c("spe filename", "Dimension Reduction", "name for output", "prelim cluster to use")
    ),
    ncol = 5
)
opt <- getopt(spec)

## get opts
dimred <- opt$dimred
name <- opt$name

message("Load spe from:", opt$spe)
message("Dimension Reduction:", dimred)
message("Output name:", name)
message("Output name:", opt$prelim)

spe_hdf5_path <- here("processed-data", "spe_objects", opt$spe)
stopifnot(file.exists(spe_hdf5_path))

## Choose k
k <- as.numeric(
    #   Only one of these environment variables will be defined, so grab the
    #   defined one (handle SGE or SLURM)
    paste0(Sys.getenv("SLURM_ARRAY_TASK_ID"), Sys.getenv("SGE_TASK_ID"))
)
k_nice <- sprintf("%02d", k)

message("k:", k_nice)

## Load the data
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here(spe_hdf5_path))
stopifnot(dimred %in% reducedDimNames(spe))

## choose prelim clusters

prelim_cluster <- NULL
if(!is.null(opt$prelim)){
    
    prelim_cluster <- sprintf("10x_kmeans_%d_clusters", k)
    message("prelim cluster col:", prelim_cluster)
    stopifnot(prelim_cluster %in% colnames(colData(spe)))

    if(opt$prelim == "10x_id"){
        spe[[prelim_cluster]] <- paste0(spe$sample_id, "_", spe[[prelim_cluster]])        
    }
    
    table(spe[[prelim_cluster]])
}


## Set Seed
set.seed(20240905)

## Create output directories
dir_plots <- here("plots", "05_spe_correct_cluster", "05_BayesSpace", paste0("clusters_BayesSpace-", name), paste0(name, "_k", k_nice))
if(!dir.exists(dir_plots)) dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

dir_rdata <- here("processed-data", "05_spe_correct_cluster", "05_BayesSpace", paste0("clusters_BayesSpace-", name))
if(!dir.exists(dir_rdata)) dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)


## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message(Sys.time(), sprintf(" - Running spatialCluster: k = %d, dimred = %s, inint = %s", k, dimred, prelim_cluster))
spe <- BayesSpace::spatialCluster(spe, 
                                  use.dimred = dimred,
                                  q = k, 
                                  platform = "Visium",
                                  init = prelim_cluster,
                                  nrep = 10000)

message(Sys.time(), " - Format and Export")

table(spe$spatial.cluster)

spe$bayesSpace_temp <- as.factor(spe$spatial.cluster)
bayesSpace_name <- paste0("BayesSpace_",name,"_k", k_nice)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here(dir_rdata, "clusters_BayesSpace")
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
    point_size = 1
)

#  Vis_clus for each sample
walk(sample_ids, function(samp){
    spot_plot <- vis_clus(
        spe = spe,
        sampleid = samp,
        clustervar = bayesSpace_name,
        # sort_clust = FALSE,
        colors = cols,
        point_size = 2
    )
    ggsave(spot_plot, filename = here(dir_plots, paste0(bayesSpace_name, "-", samp, ".pdf")))
})

## in the future, also plot reduced dims by cluster

# slurmjobs::job_single('05_BayesSpace', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
