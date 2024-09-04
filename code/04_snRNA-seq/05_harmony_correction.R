library("SingleCellExperiment")
library("harmony")
library("scater")
library("HDF5Array")
library("here")
library("sessioninfo")


## Load "uncorrected" sce data (postQC, GLMPCA)
message(Sys.time(), " - load data")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_uncorrected"))
sce

## Choose Correction
correction = "sample_id"
stopifnot(correction %in% colnames(colData(sce)))
message("Correcting by: ", correction)

## Run harmony
## needs PCA
reducedDim(sce, "PCA") <- reducedDim(sce, "GLMPCA_approx")

message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce, group.by.vars = correction, verbose = TRUE)

## Remove redundant PCA
reducedDim(sce, "PCA") <- NULL

#### TSNE & UMAP ####
set.seed(602)
message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = "HARMONY")

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY")

message("Done TSNE + UMAP - Saving data...", Sys.time())

# save(sce, file = here("processed-data", "spe_objects", paste0("sce_harmony_", correction, ".Rdata")))
save(sce, file = here("processed-data", "spe_objects", "sce_harmony.Rdata"))

# saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "03_build_sce", paste0("sce_harmony_", correction)))
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "spe_objects", "sce_harmony"))

# slurmjobs::job_single('05_harmony_correction', create_shell = TRUE, memory = '100G', command = "Rscript 05_harmony_correction.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
