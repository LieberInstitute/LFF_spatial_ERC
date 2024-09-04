library("SingleCellExperiment")
library("harmony")
library("scater")
library("HDF5Array")
library("here")
library("sessioninfo")

## Choose Correction
correction = ""

## Load empty-free sce data
message(Sys.time(), " - load data")
sce_uncorrected <- loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "sce_uncorrected"))

stopifnot(correction %in% colnames(colData(sce_uncorrected)))
message("Correcting by: ", correction)

## Run harmony
## needs PCA
reducedDim(sce_uncorrected, "PCA") <- reducedDim(sce_uncorrected, "GLMPCA_approx")

message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce_uncorrected, group.by.vars = correction, verbose = TRUE)

## Remove redundant PCA
reducedDim(sce, "PCA") <- NULL

#### TSNE & UMAP ####
set.seed(602)
message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = "HARMONY")

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY")

message("Done TSNE + UMAP - Saving data...", Sys.time())

save(sce, file = here("processed-data", "03_build_sce", paste0("sce_harmony_", correction, ".Rdata")))

# saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "03_build_sce", paste0("sce_harmony_", correction)))
saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "spe_objects", "sce_harmony"))


# slurm_jobs::job_single('05_harmony_correction', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 05_harmony_correction.R Sample")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
