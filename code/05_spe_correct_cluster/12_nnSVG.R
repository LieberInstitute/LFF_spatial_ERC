## October 2024, Louise Huuki-Myers
## Run nnSVG by sample 
## Adapted from https://github.com/LieberInstitute/visiumStitched_brain/blob/9d2f716bd76a3359b2aeca1cac3c90720c12bfb1/code/03_stitching/04_nnSVG.R#L7

# library(tidyverse)
library("SpatialExperiment")
library("HDF5Array")
library("Matrix")
library("nnSVG")
library("getopt")
library("here")
library("sessioninfo")

#### Import command-line parameters ####
spec <- matrix(
    c(
        c("sample"),
        c("s"),
        rep("1", 1),
        rep("character", 1),
        rep("sample in spe", 1)
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

message(sprintf("Sample [%s]: %s", Sys.getenv("SLURM_ARRAY_TASK_ID"), opt$sample))

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "12_nnSVG")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

##sub set to sample
spe <- spe[,spe$sample_id == opt$sample]

####   Bring assays into memory for speed ####
message(Sys.time(), " - Bringing into memory")
assays(spe) <- list(
    counts = as(assays(spe)$counts, "dgCMatrix"),
    logcounts = as(assays(spe)$logcounts, "dgCMatrix")
)

#### Filter lowly expressed and mitochondrial genes ####
set.seed(0)

message(Sys.time(), " - Filtering genes and spots")
spe <- filter_genes(
    spe,
    filter_genes_ncounts = 3,
    filter_genes_pcspots = 0.5,
    filter_mito = TRUE
)

message("Post filter spots: ", ncol(spe))
message("Post filter genes: ", nrow(spe))


#####   Run nnSVG and export results ####
message(Sys.time(), " - Running nnSVG")
spe <- nnSVG(spe)

message(Sys.time(), " - Exporting results")
write.csv(as_tibble(rowData(spe)), here(data_dir, sprintf("nnSVG_%s.csv", opt$sample)))

# slurmjobs::job_loop(loops = list(sample = sort(unique(spe$sample_id))), '12_nnSVG', create_shell = TRUE, memory = '50G')

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
