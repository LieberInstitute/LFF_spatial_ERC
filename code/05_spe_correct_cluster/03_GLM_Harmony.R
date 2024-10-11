# Aug, 2024 - Louise Huuki-Myers
# Compute GLM PCA & Run Harmony batch correction for Visium data

library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scran")
library("scater")
library("harmony")
library("BiocParallel")
library("purrr")
library("scry")
library("HDF5Array")

plot_dir <- here("plots", "05_spe_correct_cluster", "03_GLM_Harmony")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "03_GLM_Harmony")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


#### Load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))


####  Compute GLM-PCA ####
message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial", sorted = FALSE, batch = as.factor(spe$sample_id))
# spe <- devianceFeatureSelection(spe, assay = "counts", fam = "poisson", sorted = FALSE, batch = as.factor(spe$sample_id)) 

## plot binomial and poison deviance in first 100 selected genes ?

## calculate residuals from binomial model

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals( # default params
    spe,
    assay = "counts",
    fam = "binomial",
    type = "deviance"
    # batch = as.factor(spe$sample_id) # batch membership of observations - use "round" ?
)

## Get HDGs
top.hdgs <- map(c(`1k` = 1000, `2k` = 2000, `5k` = 5000),
                ~rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:.x]
)

save(top.hdgs, file = here(data_dir, "top_hdgs.rdata"))

walk2(top.hdgs, names(top.hdgs), function(hdgs, name_hdgs){
    message(Sys.time(), " - Running GLM-PCA: ", name_hdgs)
    spe <<- runPCA(
        spe,
        exprs_values = "binomial_deviance_residuals",
        subset_row = hdgs,
        ncomponents = 100,
        name = paste0("GLMPCA_approx_", name_hdgs),
        BSPARAM = BiocSingular::IrlbaParam()
    )
    
})


#### Harmony Batch Correction ####

## Harmony uses the PCA slot
reducedDim(spe, "PCA") <- reducedDim(spe, "GLMPCA_approx_2k")

message(Sys.time(), " - Harmony")
spe <- RunHarmony(spe, "sample_id")

message(Sys.time(), " - Harmony UMAP")
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")

colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

message(Sys.time(), " - HARMONY TSNE")
spe <- runTSNE(spe,
               dimred = "HARMONY",
               name = "TSNE.HARMONY")

## Remove redundant PCA
reducedDim(spe, "PCA") <- NULL

message(Sys.time(), " - Saving HDF5 SPE")
saveHDF5SummarizedExperiment(
    spe,
    dir = here(data_dir, "spe_GLM_Harmony"), replace = TRUE
)

# slurmjobs::job_single('03_GLM_Harmony', create_shell = TRUE, memory = '25G', command = "Rscript 03_GLM_Harmony.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()