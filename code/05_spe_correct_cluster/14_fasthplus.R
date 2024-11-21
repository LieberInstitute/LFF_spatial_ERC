## October 2024, Louise Huuki-Myers
## Evaluate fast h+ metric for spatial clusters
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/06_fasthplus/01_fasthplus.R

library("fasthplus")
library("spatialLIBD")
library("here")
library("sessioninfo")
library("getopt")

data_dir <- here("processed-data", "05_spe_correct_cluster", "14_fasthplus")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Import command-line parameters
spec <- matrix(
    c(  "cluster", "i", "1", "character", "Name of cluster",
        "k", "k", "1", "numeric", "Number of clusters",
        "drop_WM", "d", "1", "boolean", "drop White Matter - defined in script"
        ),
    ncol = 5, byrow = TRUE
)
opt <- getopt(spec)

# for testing
# opt <- list(cluster ="BayesSpace_SVGm", k = 10, drop_WM = TRUE)

print("Using the following parameters:")
print(opt)

nice_k = sprintf("k%02d",opt$k)
cluster_k <- sprintf("%s_k%02d", opt$cluster ,opt$k)
message("Runnign fast h+ on: ", cluster_k)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))
stopifnot(cluster_k %in% colnames(colData(spe)))

# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
dim(reducedDims(spe)$HARMONY)
# [1] 122202     50

## drop WM
if(opt$drop_WM){
    message("Dropping WM")
    spe <- spe[, spe$BayesSpace_SVGm_k02 != "Sp02D02"]
} 

#### Set up eval ####
message(Sys.time(), " - Setup fasthplus")

## find intitial t
find_t <- function(L, proportion = 0.05) {
    initial_t <- floor(length(L) * proportion)
    smallest_cluster_size <- min(table(L))
    n_labels <- length(unique(L))
    ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

t <- find_t(L = colData(spe)[[cluster_k]], proportion = 0.01)
message("Initial t: ", t)

message("cluster proportions:")
(cluster_prop <- table(colData(spe)[[cluster_k]])/ncol(spe))

min_prop <- 0.01/opt$k

message("Testing for clusters < ", min_prop)
bad_clusters <- which(cluster_prop < min_prop)

if (length(bad_clusters) > 0) {
    message("For ", nice_k, " drop small clusters: ", paste(names(bad_clusters), collapse = ", "))
    spe <- spe[, !colData(spe)[[cluster_k]] %in% names(bad_clusters)]
    t <- find_t(colData(spe)[[cluster_k]], 0.01)
    message("Updated t: ", t)
} else message("No small clusters")


#### Run fasthplus ####
set.seed(202310)

message(Sys.time(), " - Run fasthplus")
fasthplus <- hpb(D = reducedDims(spe)$HARMONY, 
                 L = colData(spe)[[cluster_k]],
                 t = t, 
                 r = 30)

message(Sys.time(), " - Done")
message("Results: ", cluster_k,  "H+ :", fasthplus)

## save results
results <- data.frame(cluster = cluster_k, fasthplus = fasthplus)
write.table(results, file = here(data_dir, sprintf("fasthplus_results-%s.csv", opt$cluster), append = TRUE))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()