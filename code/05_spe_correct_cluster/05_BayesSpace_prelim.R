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

## to test
# opt <- list(spe = "spe_postQC", dimred = "HARMONY", name = "PCA_Harmony_prelim", prelim = "10x")

## get opts
dimred <- opt$dimred
name <- opt$name

message("Load spe from:", opt$spe)
message("Dimension Reduction:", dimred)
message("Prelim input:", opt$prelim)
message("Output name:", name)

spe_hdf5_path <- here("processed-data", "spe_objects", opt$spe)
stopifnot(file.exists(spe_hdf5_path))

## Choose k
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
k_nice <- sprintf("%02d", k)

message("k:", k)

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

## Create output directories
dir_plots <- here("plots", "05_spe_correct_cluster", "05_BayesSpace", paste0("clusters_BayesSpace-", name), paste0(name, "_k", k_nice))
if(!dir.exists(dir_plots)) dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

dir_rdata <- here("processed-data", "05_spe_correct_cluster", "05_BayesSpace", paste0("clusters_BayesSpace-", name))
if(!dir.exists(dir_rdata)) dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

#### Cluster ####
message(Sys.time(), " - Set Up Clustering")

## Set Seed
set.seed(20240905)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

## Fix colnames for row and col
spe$row <- spe$array_row
spe$col <- spe$array_col

type(spe[[prelim_cluster]]) #[1] "integer"
spe$prelim <- as.double(spe[[prelim_cluster]])
type(spe[["prelim"]]) 

message(Sys.time(), sprintf(" - Running spatialCluster: k = %d, dimred = %s, inint = %s", k, dimred, prelim_cluster))
spe <- BayesSpace::spatialCluster(spe, 
                                  use.dimred = dimred,
                                  q = k, 
                                  # init = prelim_cluster,
                                  init = spe$prelim,
                                  nrep = 10000)

# Error: Not compatible with requested type: [type=character; target=double]. # "init" as char, double, or int 
# Traceback 
# 4: stop(structure(list(message = "Not compatible with requested type: [type=character; target=double].", 
#                        call = NULL, cppstack = NULL), class = c("Rcpp::not_compatible", 
#                                                                 "C++Error", "error", "condition")))
# 3: cluster.FUN(Y = as.matrix(Y), df_j = df_j, nrep = nrep, n = n, 
#                d = d, gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, 
#                alpha = alpha, beta = beta)
# 2: cluster(Y, q, df_j, init = init, model = model, precision = precision, 
#            mu0 = mu0, lambda0 = lambda0, gamma = gamma, alpha = alpha, 
#            beta = beta, nrep = nrep)
# 1: BayesSpace::spatialCluster(spe, use.dimred = dimred, q = k, init = "prelim", 
#                               nrep = 10000)

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
# walk(sample_ids, function(samp){
#     spot_plot <- vis_clus(
#         spe = spe,
#         sampleid = samp,
#         clustervar = bayesSpace_name,
#         # sort_clust = FALSE,
#         colors = cols,
#         point_size = 2
#     )
#     ggsave(spot_plot, filename = here(dir_plots, paste0(bayesSpace_name, "-", samp, ".pdf")))
# })

## in the future, also plot reduced dims by cluster

# slurmjobs::job_single('05_BayesSpace', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
