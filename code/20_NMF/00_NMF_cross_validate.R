## Louise Huuki-Myers, Jan 2026
## Find rank that minimizes mean squared error for the data

#### Set up ####
library("singlet")
library("SingleCellExperiment")
library("RcppML")
library("sessioninfo")
library("here")

data_dir <- here("processed-data", "20_NMF", "00_NMF_cross_validate")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "00_NMF_cross_validate")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

assay(sce,'logcounts')

## Drop Br1289
sce <- sce[, sce$BrNum != "Br1289"]

## logcounts to sparse matrix
message(Sys.time(), " - logcounts to sparse matrix")
logcounts(sce) <- as(logcounts(sce), "sparseMatrix")

## Run cross_validate_nmf
message(Sys.time(), " - Run cross_validate_nmf")

cvnmf <- cross_validate_nmf(
    logcounts(sce),   
    ranks=c(5,10,20, 30, 40, 50, 75, 100, 125),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)


# Save NMF results
message(Sys.time(), " - Done NMF...Saving")
saveRDS(cvnmf, file = here(data_dir, "ERC_cvnmf.RDS"))

# Plot cross validation
pdf(here(plot_dir,"nmf_cv_results_lower_k.pdf"))

plot(cvnmf)

dev.off()


# slurmjobs::job_single('00_NMF_cross_validate', create_shell = TRUE, memory = '200G', command = "Rscript 00_NMF_cross_validate.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()