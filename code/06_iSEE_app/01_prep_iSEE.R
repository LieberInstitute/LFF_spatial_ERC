library("SingleCellExperiment")
library("here")
library("lobstr")
library("sessioninfo")
library("HDF5Array")
library("Matrix")
library("tidyverse")

## Load HD5F sce
message(Sys.time(), "- load SCEe")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "sce_objects", "sce_ERC_subcluster"))
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

# # sourcing official color palette
sn_colors <- metadata(sce)$cell_type_colors

## Drop metadata we don't need
metadata(sce) <- list()
sce$path <- NULL
sce$total <- NULL

## Check final size
lobstr::obj_size(sce)
# 5.47 GB

rownames(sce) <- rowData(sce)$gene_name

#### Add MeanRatio Marker Gene Details ####
load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio", "marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)
marker_stats_MeanRatio |> dplyr::count(cellType.target)

marker_anno <- marker_stats_MeanRatio |>
    filter(MeanRatio.rank <= 50 & MeanRatio > 1) |>
    select(gene,
           cellType.target,
           MeanRatio.rank,
           MeanRatio,
           MeanRatio.anno) |>
    column_to_rownames("gene")

any(duplicated(rownames(marker_anno)))

rowData(sce) <- cbind(rowData(sce), marker_anno[rownames(sce),])

rowData(sce)[which(rowData(sce)$MeanRatio.rank ==1),]


sce
# class: SingleCellExperiment 
# dim: 38606 140119 
# metadata(0):
#     assays(1): logcounts
# rownames(38606): ENSG00000290825 ENSG00000243485 ... ENSG00000278817 ENSG00000277196
# rowData names(4): ID Symbol Type binomial_deviance
# colnames(140119): cell1 cell2 ... cell160742 cell160743
# colData names(34): Barcode sample_id ... low_detected_batch sizeFactor

## save RDS
message(Sys.time(), "- Save data RDS")
saveRDS(sce, file = here("code", "06_iSEE_app", "sce_ERC_iSEE.rds"))
saveRDS(sn_colors, file = here("code", "06_iSEE_app", "sn_colors.rds"))

## save HDF5
message(Sys.time(), "- Save data HDF5")
saveHDF5SummarizedExperiment(sce, dir = here("code", "06_iSEE_app", "sce_ERC_iSEE"), replace=TRUE)


# slurmjobs::job_single('01_prep_iSEE', create_shell = TRUE, memory = '25G', command = "Rscript 01_prep_iSEE.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
