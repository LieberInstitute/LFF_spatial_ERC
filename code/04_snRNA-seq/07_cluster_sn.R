## Louise Huuki-Myers, September 2024
## use single nearest neighbors + walk trap to cluster single nuc data

library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("HDF5Array")
library("here")
library("sessioninfo")

## Prep directories
data_dir <- here("processed-data", "04_snRNA-seq", "07_cluster_sn")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "sce_harmony"))
sce

## define k
k=20

## set seed
set.seed(20240905)

## Build SNN graph
message(Sys.time(), " - running buildSNNGraph: k =", k)
snn.gr <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")

## Run walk trap clustering
message(Sys.time(), " - running walktrap")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters)

## save data
message(Sys.time() , " - saving data")
save(clusters, file = here(data_dir, sprintf("walktrap_snn_k%02d_clusters.Rdata", k)))

## Save final sce w/ annotations
# save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# slurmjobs::job_single('07_cluster_sn', create_shell = TRUE, memory = '100G', command = "Rscript 07_cluster_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
