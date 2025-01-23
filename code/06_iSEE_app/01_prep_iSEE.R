library("SingleCellExperiment")
library("here")
library("lobstr")
library("sessioninfo")
library("HDF5Array")
library("Matrix")

## Load HD5F sce
message(Sys.time(), "- load Harmony corrected sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))
sce

colnames(colData(sce))

## Drop data we don't need for iSEE
assayNames(sce)
# [1] "counts"                      "binomial_deviance_residuals" "logcounts"
assays(sce)$counts <- NULL
assays(sce)$binomial_deviance_residuals <- NULL

class(logcounts(sce))
# [1] "DelayedMatrix"

message(Sys.time(), "- Convert logcounts to sparse Matrix")
logcounts(sce) <- as(logcounts(sce), "sparseMatrix")  

## Drop metadata we don't need
metadata(sce) <- list()
sce$path <- NULL
sce$total <- NULL

# # sourcing official color palette
# source(
#     file = here("code", "99_paper_figs", "source_colors.R"),
#     echo = TRUE
# )

## Check final size
lobstr::obj_size(sce)
# 5.87 GB

sce
# class: SingleCellExperiment 
# dim: 38606 140119 
# metadata(0):
#     assays(1): logcounts
# rownames(38606): ENSG00000290825 ENSG00000243485 ... ENSG00000278817 ENSG00000277196
# rowData names(4): ID Symbol Type binomial_deviance
# colnames(140119): cell1 cell2 ... cell160742 cell160743
# colData names(34): Barcode sample_id ... low_detected_batch sizeFactor

sce 
saveRDS(sce, file = here("code", "06_iSEE_app", "sce_ERC_iSEE.rds"))
# saveRDS(sn_colors, file = here("code", "06_iSEE_app", "sn_colors.rds"))

# slurmjobs::job_single('01_prep_iSEE', create_shell = TRUE, memory = '25G', command = "Rscript 01_prep_iSEE.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
