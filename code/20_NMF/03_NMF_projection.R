## Louise Huuki-Myers, Jan 2026
## Correlate NMF factors with variables

#### Set up ####

# library("SpatialExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("ComplexHeatmap")
library("spatialLIBD")
library("projectR")

data_dir <- here("processed-data", "20_NMF", "03_NMF_projection")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "03_NMF_projection")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####

## snRNA-seq
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## spatail data
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

## Drop Br1289
sce <- sce[, sce$BrNum != "Br1289"]
spe <- spe[, sce$BrNum != "Br1289"]

#### NMF data ####
nmf <- readRDS(here("processed-data", "20_NMF", "01_runNMF", "ERC_sn_nmf.RDS"))

# extract patterns
patterns <- t(nmf@h)

loadings <- nmf@w
rownames(loadings) <- rownames(sce)
loadings[1:5,1:5]

## project loadings to spatial data
i <- intersect(rowData(spe)$gene_id,rownames(loadings))
loadings <- loadings[rownames(loadings) %in% i,]
spe <- spe[rowData(spe)$gene_id %in% i,]
loadings <- loadings[match(rowData(spe)$gene_id,rownames(loadings)),]

logcounts <- logcounts(spe)
#loadings <- as(loadings, "dgCMatrix")

# proj <- project(w=loadings, data=logcounts)

proj <- projectR(data = as(logcounts, "Matrix"), loadings = loadings)

# remove rowSums == 0
proj <- proj[rowSums(proj) != 0,]

proj <- apply(proj,1,function(x){x/sum(x)})

saveRDS(proj, file = here(data_dir, "ECR_nmf_projection.Rds"))

# slurmjobs::job_single('03_NMF_projection', create_shell = TRUE, memory = '200G', command = "Rscript 04_multistudy_Oligo_subcluster.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
