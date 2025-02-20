# Aug, 2024 - Louise Huuki-Myers
# Read in postQC SCE created in 02_droplet_qc
# Run GLM PCA

library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("HDF5Array")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "04_snRNA-seq", "03_GLM_PCA")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## Load HD5F
message(Sys.time(), " - load post-QC sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "hdf5_sce_postQC"))
sce

#### Deviance Feat Selection ####
set.seed(606)

message(Sys.time(), " - running Deviance Feat. Selection")
sce <- devianceFeatureSelection(sce,
    assay = "counts", fam = "binomial", sorted = F,
    batch = as.factor(sce$seq_round)
)

# This temp file just used for getting batch-corrected components (drops a variety of entries)

pdf(here(plot_dir, "binomial_deviance.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
    type = "l", xlab = "ranked genes",
    ylab = "binomial deviance", main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()

message(Sys.time(), " - running nullResiduals - ")
sce <- nullResiduals(sce,
    assay = "counts", 
    fam = "binomial", # default params
    type = "deviance"
)


hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]

message(Sys.time(), " - running PCA - ")
sce <- runPCA(sce,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs,
    ncomponents = 100,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

#### Reduce Dims####

message(Sys.time(), " - running TSNE")
sce <- runTSNE(sce, dimred = "GLMPCA_approx")

message(Sys.time(), " - running UMAP")
sce <- runUMAP(sce, dimred = "GLMPCA_approx")


message("SCE contents:")
sce

message("Reduced Dim Names:")
reducedDimNames(sce)

data_dir <- here("processed-data", "sce_objects", "sce_uncorrected")
message(Sys.time(), " - Saving Data to: ", data_dir)
saveHDF5SummarizedExperiment(sce, dir = data_dir, replace = TRUE)

# slurmjobs::job_single('03_GLM_PCA', create_shell = TRUE, memory = '100G', command = "Rscript 03_GLM_PCA.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
