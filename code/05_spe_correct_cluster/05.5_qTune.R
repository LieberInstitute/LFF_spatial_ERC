## Louise Huuki-Myers, Nov 2024
## Run BayesSpace qTune on selected reduced dim for k2-28

## Required libraries
library("spatialLIBD")
library("BayesSpace")
library("Polychrome")
library("tidyverse")
library("HDF5Array")
library("getopt")
library("scater")
library("harmony")
library("here")
library("sessioninfo")

## get input file from params
spec <- matrix(
    c(
        c("spe", "dimred", "name"),
        c("s", "d", "n"),
        c("1", "2", "3"),
        c("character", "character", "character"),
        c("spe filename", "Dimension Reduction", "name for output")
    ),
    ncol = 5
)
opt <- getopt(spec)

## get opts
dimred <- opt$dimred
name <- opt$name

message("Load spe from:", opt$spe)
message("Dimension Reduction:", dimred)
message("Output name:", name)

## Create output directories
plot_dir <- here("plots", "05_spe_correct_cluster", "05.5_qTune")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "05.5_qTune")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")

spe_hdf5_path <- here("processed-data", "spe_objects", opt$spe)
stopifnot(file.exists(spe_hdf5_path))

#### Cluster Prep ####
message(Sys.time(), " - Set Up Clustering")
## Set Seed
set.seed(20240905)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

## Fix colaname for row and col
spe$row <- spe$array_row
spe$col <- spe$array_col

## pick number of clusters
## qTune uses PCA - replace w/ HARMONY reduced dims
message("Use ", dimred, "as 'PCA'")
reducedDim(spe, "PCA") <- reducedDim(spe, dimred) 

message(Sys.time(), " - qTune")
spe <- qTune(spe, qs=seq(2, 28), platform="Visium")

attr(spe, "q.logliks")

## save model's log likelihood
message(Sys.time(), " - Done qTune, save and plot")
write.csv(attr(spe, "q.logliks"), file = here(data_dir, sprintf("qTune_logliks_%s.csv", name)))

# qplot
pdf(here(plot_dir, sprintf("BayesSpace_erc_qplot_%s.pdf", name)))
qPlot(spe)
dev.off()

# slurmjobs::job_single('05.5_qTune', create_shell = TRUE, memory = '25G', command = "Rscript 05.5_qTune.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
