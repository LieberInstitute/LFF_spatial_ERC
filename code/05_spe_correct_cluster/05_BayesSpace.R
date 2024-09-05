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

set.seed(20240905)

## Choose k
k <- as.numeric(
    #   Only one of these environment variables will be defined, so grab the
    #   defined one (handle SGE or SLURM)
    paste0(Sys.getenv("SLURM_ARRAY_TASK_ID"), Sys.getenv("SGE_TASK_ID"))
)
k_nice <- sprintf("%02d", k)

## Create output directories
dir_plots <- here("plots", "05_spe_correct_cluster", "05_BayesSpace", paste0("k", k_nice))
if(!dir.exists(dir_plots)) dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

dir_rdata <- here("processed-data", "05_spe_correct_cluster", "05_BayesSpace")
if(!dir.exists(dir_rdata)) dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
if(!dir.exists(here(dir_rdata, "clusters_BayesSpace"))) dir.create(here(dir_rdata, "clusters_BayesSpace"), showWarnings = FALSE, recursive = TRUE)

## Load the data
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_GLM_Harmony"))

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message(Sys.time(), " - Running spatialCluster(): k=", k)
spe <- BayesSpace::spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)

message(Sys.time(), " - Format and Export")

table(spe$spatial.cluster)

spe$bayesSpace_temp <- as.factor(spe$spatial.cluster)
bayesSpace_name <- paste0("BayesSpace_harmony_k", k_nice)
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

# slurmjobs::job_single('05_BayesSpace', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
