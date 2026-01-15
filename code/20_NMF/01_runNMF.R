## Louise Huuki-Myers, Jan 2026
## Run NMF on ERC snRNA-seq data
## following https://github.com/LieberInstitute/spatialdACC/blob/main/code/snRNA-seq/06_NMF/01_NMF_factor.R

#### Set up ####

library("singlet")
library("dplyr")
library("ggplot2")
library("sessioninfo")
library("here")
library("SingleCellExperiment")
library("RcppML")

data_dir <- here("processed-data", "20_NMF", "01_runNMF")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("plots", "20_NMF", "01_runNMF")
# if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

assay(sce,'logcounts')[1:5, 1:5]

## Run NMF
message(Sys.time(), " - Run NMF")
options(RcppML.threads=4)
x <- RcppML::nmf(assay(sce,'logcounts'),
                 k=75,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 114,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE
)

# Save NMF results
message(Sys.time(), " - Done NMF...Saving")
saveRDS(x, file = here(data_dir, "ERC_sn_nmf.RDS"))

# slurmjobs::job_single('01_runNMF', create_shell = TRUE, memory = '100G', command = "Rscript 01_runNMF")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


