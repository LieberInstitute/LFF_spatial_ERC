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
message(Sys.time(), "- load post-QC sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "hdf5_sce_postQC"), prefix = "")
dim(sce)

#### DEviance Feat Selection ####
set.seed(606)

message("running Deviance Feat. Selection - ", Sys.time())
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

message("running nullResiduals - ", Sys.time())
sce <- nullResiduals(sce,
    assay = "counts", fam = "binomial", # default params
    type = "deviance"
)


hdgs.hb <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]

message("running PCA - ", Sys.time())
sce_uncorrected <- runPCA(sce,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb, ncomponents = 100,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

#### Reduce Dims####

message("running TSNE - ", Sys.time())
sce_uncorrected <- runTSNE(sce_uncorrected, dimred = "GLMPCA_approx")

message("running UMAP - ", Sys.time())
sce_uncorrected <- runUMAP(sce_uncorrected, dimred = "GLMPCA_approx")

message("Saving Data - ", Sys.time())
saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "03_build_sce", "sce_uncorrected"))

# slurmjobs::job_single('03_GLM_PCA', create_shell = TRUE, memory = '100G', command = "Rscript 03_GLM_PCA.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
