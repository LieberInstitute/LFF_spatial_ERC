## Louise Huuki-Myers, September 2024
## use single nearest neighbors + walk trap to find prelim cluster single nuc data

library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("HDF5Array")
library("getopt")
library("here")
library("sessioninfo")

# Import command-line parameters
spec <- matrix(
    c(
        c("k"),
        c("k"),
        rep("1", 1),
        rep("integer", 1),
        rep("k for buildSNNGraph", 1)
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

k = opt$k

## Prep directories
data_dir <- here("processed-data", "04_snRNA-seq", "07_cluster_sn_prelim")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "sce_harmony"))
sce


## set seed
set.seed(20240905)

## Build SNN graph
message(Sys.time(), " - running buildSNNGraph: k =", k)
snn.gr <- buildSNNGraph(sce, k = k, use.dimred = "HARMONY")

## Run walk trap clustering
message(Sys.time(), " - running walktrap")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters)

## save data
message(Sys.time() , " - saving data")
save(clusters, file = here(data_dir, sprintf("walktrap_snn_k%02d_clusters_prelim.Rdata", k)))

# slurmjobs::job_single('07_cluster_sn_prelim', create_shell = TRUE, memory = '100G', command = "Rscript 07_cluster_sn_prelim.R -k 10")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
