
library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(HDF5Array)
library(Matrix)
library(nnSVG)

spe_dir <- here("processed-data", "03_stitching", "spe")
out_path <- here(
    "processed-data", "03_stitching", "nnSVG_out", "Br2719.csv"
)

set.seed(0)
dir.create(dirname(out_path), showWarnings = FALSE)

message(Sys.time(), " | Loading SpatialExperiment")
spe <- loadHDF5SummarizedExperiment(spe_dir)

#-------------------------------------------------------------------------------
#   Bring into memory for speed
#-------------------------------------------------------------------------------

message(Sys.time(), " | Bringing into memory")
assays(spe) <- list(
    counts = as(assays(spe)$counts, "dgCMatrix"),
    logcounts = as(assays(spe)$logcounts, "dgCMatrix")
)

#-------------------------------------------------------------------------------
#   Filter lowly expressed and mitochondrial genes
#-------------------------------------------------------------------------------

message(Sys.time(), " | Filtering genes and spots")
spe <- filter_genes(
    spe,
    filter_genes_ncounts = 3,
    filter_genes_pcspots = 0.5,
    filter_mito = TRUE
)
message("Dimensions of spe after filtering:")
print(dim(spe))

#-------------------------------------------------------------------------------
#   Run nnSVG and export results
#-------------------------------------------------------------------------------

message(Sys.time(), " | Running nnSVG")
spe <- nnSVG(spe)

message(Sys.time(), " | Exporting results")
write_csv(as_tibble(rowData(spe)), out_path)

session_info()
