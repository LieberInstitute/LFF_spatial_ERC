## Louise Huuki-Myers, February 2025
## re-process data after luster QC
## drop low quality prelim-clusters, re-run normalization, feature selection, dimension reduction

## Required libraries
library("SingleCellExperiment")
library("tidyverse")
library("HDF5Array")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")

## source functions
source(here("code", "utils", "my_plot_reduced_dim_ALL.R"))

## Set up dirs
data_dir <- here("processed-data", "04_snRNA-seq", "10_reprocess_quality_sn")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "10_reprocess_quality_sn")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

## load cluster outputs
load(here("processed-data", "04_snRNA-seq", "07_cluster_sn_prelim", "walktrap_snn_k10_clusters_prelim.Rdata"), verbose = TRUE)
sce$snn_k10 <- sprintf("k10c%02d", clusters)

## load cluster QC output
cluster_qc <- read.csv(here("processed-data", "04_snRNA-seq", "09_cluster_QC", "sn_clusterQC_info.csv"))
clusters_passQC <- cluster_qc$snn_k10[cluster_qc$pass_clusterQC]

#### Filter SCE ####
message("Nuclei pre-cluster QC: ", ncol(sce))
sce <- sce[, sce$snn_k10 %in% clusters_passQC]
message("Nuclei post-cluster QC: ", ncol(sce))

#### Normalize logcounts ####
message(Sys.time(), " - Quick Cluster")
sce$quick_cluster <- quickCluster(sce) 

message(Sys.time(), " - logNormCounts")
sce <- computeSumFactors(sce, cluster=sce$quick_cluster)
sce <- logNormCounts(sce)

#### GLM PCA ####
set.seed(219)

message(Sys.time(), " - running Deviance Feat. Selection")
sce <- scry::devianceFeatureSelection(sce,
                                assay = "counts", fam = "binomial", sorted = F,
                                batch = as.factor(sce$seq_round)
)

## plot biononial deviance distribution
pdf(here(plot_dir, "sn_reprocess_binomial_deviance.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance", main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()

message(Sys.time(), " - running nullResiduals")
sce <- nullResiduals(sce,
                     assay = "counts", 
                     fam = "binomial", # default params
                     type = "deviance"
)


hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]

message(Sys.time(), " - running PCA")
sce <- runPCA(sce,
              exprs_values = "binomial_deviance_residuals",
              subset_row = hdgs,
              ncomponents = 100,
              name = "GLMPCA_approx",
              BSPARAM = BiocSingular::IrlbaParam()
)

#### Reduce Dims - UNCORRECTED ####

message(Sys.time(), " - running TSNE")
sce <- runTSNE(sce, dimred = "GLMPCA_approx")

message(Sys.time(), " - running UMAP")
sce <- runUMAP(sce, dimred = "GLMPCA_approx")

message("SCE contents:")
sce

message("Reduced Dim Names:")
reducedDimNames(sce)

#### plot uncorrected reduced dims ####
message(Sys.time(), " - Reduced Dim plots - UNCORRECTED")
my_plot_reduced_dim_ALL(prefix = "ERC_sn_reprocess", suffix = "uncorrected")

#### Batch Correction ####
correction = "sample_id"
message("HARMONY Correcting by: ", correction)

## Run harmony
## needs PCA
reducedDim(sce, "PCA") <- reducedDim(sce, "GLMPCA_approx")

message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce, group.by.vars = correction, verbose = TRUE)

## Remove redundant PCA
reducedDim(sce, "PCA") <- NULL

#### Reduce Dims - CORRECTED ####
set.seed(602)
message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = "HARMONY")

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY")

#### Save data ####
message("Done TSNE + UMAP - Saving data...", Sys.time())
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce_objects", "sce_reprocess"), replace = TRUE)

#### plot corrected reduced dims ####
message(Sys.time(), " - Reduced Dim plots - CORRECTED")
my_plot_reduced_dim_ALL(prefix = "ERC_sn_reprocess", suffix = "uncorrected")

# slurmjobs::job_single('10_reprocess_quality_sn', create_shell = TRUE, memory = '100G', command = "Rscript 10_reprocess_quality_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

