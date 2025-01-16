## Louise Huuki-Myers, January 2025
## Preform hierarchical clustering on ssn clusters

## load libraries
library("SingleCellExperiment")
library("HDF5Array")
library("scater")
library("jaffelab")
library("dendextend")
library("dynamicTreeCut")
library("here")
library("sessioninfo")
library("dplyr")
library("getopt")

# Import command-line parameters
spec <- matrix(
    c(
        c("ssn"),
        c("s"),
        rep("1", 1),
        rep("character", 1),
        rep("ssn resolution", 1)
    ),
    ncol = 5
)
opt <- getopt(spec)

#test
# opt <- list(snn = "k10")

print("Using the following parameters:")
print(opt)

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "08_hierarchical_cluster")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "08_hierarchical_cluster")
if(!dir.exists(data_dir)) dir.create(data_dir)

## load cluster data
cluster_fn <- here("processed-data", "04_snRNA-seq", "07_cluster_sn", sprintf("walktrap_snn_%s_clusters.Rdata", opt$snn))
stopifnot(file.exists(cluster_fn))

load(cluster_fn, verbose = TRUE)
#clusters

## Load sce object
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_harmony"))
sce

# Assign as 'prelimCluster'
sce$prelimCluster <- sprintf("%sc%02d", opt$snn, clusters)

message(opt$snn, ": n clusters: ", length(unique(sce$prelimCluster)))
table(sce$prelimCluster)

clusIndexes <- splitit(sce$prelimCluster)

message(Sys.time(), " - Pseudobulk")
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii) {
    rowSums(assays(sce)$counts[, ii])
})

message("PB dim")
dim(prelimCluster.PBcounts)
# [1] 36601   297

message("Zero Sum Genes")
table(rowSums(prelimCluster.PBcounts) == 0)

message(Sys.time(), " - Get Lib Size Factors")
# Compute LSFs at this level
sizeFactors.PB.all <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
message(Sys.time() , " - Normalize", )
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {
    log2(x / sizeFactors.PB.all + 1)
}))

## Perform hierarchical clustering
message(Sys.time(), " - Cluster Again")
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)

message(Sys.time(), " - Plot")
# Print for future reference
pdf(here(plot_dir, sprintf("dend_%s.pdf", opt$snn)), height = 12)
par(cex = 0.6, font = 2)
plot(dend, main = "hierarchical cluster dend", horiz = TRUE)
# abline(v = 525, lty = 2)
dev.off()

## Save data
message(Sys.time(), " - Save")
save(dend, tree.clusCollapsed, dist.clusCollapsed, file = here(data_dir, sprintf("HC_dend_%s.Rdata", opt$snn)))

# slurmjobs::job_loop(loops = list(ssn = c("k10", "k15", "k20")),
#                     name = "10_hierarchical_cluster",
#                     create_shell = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
