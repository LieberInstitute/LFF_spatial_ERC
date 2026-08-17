## Louise Huuki-Myers, Nov 2024
## Modified for seed-sensitivity testing, 2026
##
## Run BayesSpace qTune on selected reduced dim for k2-28, repeated across
## multiple seeds via SLURM array job, to assess how much the q.logliks curve
## (and specifically the q=10 vs q=11 comparison) varies with MCMC/init
## stochasticity vs. being a stable feature of the data.
##
## Seed for array task i = 20240905 + (i - 1), so:
##   task 1 -> seed 20240905  (identical to the original BayesSpace clustering run
##                             -- lets you directly confirm this script reproduces
##                             the existing qTune_logliks_SVGm.csv before trusting
##                             the other tasks' spread)
##   task 2 -> seed 20240906
##   task 3 -> seed 20240907
##   ... etc.

## Required libraries
library("spatialLIBD")
library("BayesSpace")
library("HDF5Array")
library("getopt")
library("here")
library("sessioninfo")

## get input file from params
spec <- matrix(
    c(
        c("spe", "dimred", "name"),
        c("s", "d", "n"),
        c("append", "append", "append"),
        c("character", "character", "character"),
        c("spe filename", "Dimension Reduction", "name for output")
    ),
    ncol = 5
)
opt <- getopt(spec)

## get opts
dimred <- opt$dimred
name <- opt$name

## Determine seed from SLURM array task ID
## (defaults to task 1 / seed 20240905 if run outside of an array job)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
seed <- 20240905 + (task_id - 1)

message("Load spe from:", opt$spe)
message("Dimension Reduction:", dimred)
message("Output name:", name)
message("SLURM_ARRAY_TASK_ID: ", task_id)
message("Using seed: ", seed)

## Create output directories
plot_dir <- here("plots", "05_spe_correct_cluster", "05.5_qTune", "seedTest")
if (!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "05.5_qTune", "seedTest")
if (!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")

spe_hdf5_path <- here("processed-data", "spe_objects", opt$spe)
stopifnot(file.exists(spe_hdf5_path))

spe <- HDF5Array::loadHDF5SummarizedExperiment(here(spe_hdf5_path))
stopifnot(dimred %in% reducedDimNames(spe))

#### Cluster Prep ####
message(Sys.time(), " - Set Up Clustering")
## Set Seed -- varies by array task; task 1 preserves the original 20240905
set.seed(seed)

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

## do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)

spe$array_row <- colData(spe)$array_row + auto_offset_row[spe$sample_id] ## this sets vertical stacking


## pick number of clusters
## qTune uses PCA - replace w/ HARMONY reduced dims
message("Use '", dimred, "' as 'PCA'")
reducedDim(spe, "PCA") <- reducedDim(spe, dimred)

message(Sys.time(), " - qTune (seed ", seed, ", array task ", task_id, ")")
spe <- qTune(spe, qs = seq(2, 28), platform = "Visium")

attr(spe, "q.logliks")

## save model's log likelihood
## filename includes seed + task so parallel array tasks writing to the same
## shared output directory never clobber each other
out_csv <- here(data_dir, sprintf("qTune_logliks_%s_seed%d_task%02d.csv", name, seed, task_id))
message(Sys.time(), " - Done qTune, save and plot to ", out_csv)
write.csv(attr(spe, "q.logliks"), file = out_csv)

# qplot
pdf(here(plot_dir, sprintf("BayesSpace_erc_qplot_%s_seed%d_task%02d.pdf", name, seed, task_id)))
qPlot(spe)
dev.off()

# slurmjobs::job_single('05.5_qTune_seedTest', create_shell = TRUE, memory = '50G', command = "Rscript 05.5_qTune_seedTest.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
