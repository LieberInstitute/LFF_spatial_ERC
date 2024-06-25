
library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scran")
library("scater")
library("harmony")
library("BiocParallel")
library("purrr")
library("scry")

plot_dir <- here("plots", "05_spe_correct_cluster", "01_preprocess_Harmony")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe", "01_preprocess_Harmony")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

spe <- readRDS(here("processed-data", "02_build_spe", "spe.rds"))
dim(spe)

#### log counts ####

## quick cluster
message(Sys.time(), " - Running quickCluster()")
set.seed(20240618)
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(4),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(4)
)

table(spe$scran_quick_cluster)

## plot quick clusters?

## compute size factors
message(Sys.time(), " - Run computeSumFactors()")

spe <- computeSumFactors(spe,
                         clusters = spe$scran_quick_cluster ,
                         BPPARAM = MulticoreParam(4)
)

message("summary sizeFactors")
summary(sizeFactors(spe))

# plot deconvolution size factor for each cell compared to the equivalent size factor derived from the library size
lib.sf <- librarySizeFactors(spe)

pdf(here(plot_dir, "sizeFactor.pdf"))

hist(log10(lib.sf), xlab="Log10[Size factor]", col='grey80')

plot(lib.sf, sizeFactors(spe), xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(spe$sizeFactor)))
abline(a=0, b=1, col="red")

dev.off()

spe <- logNormCounts(spe)

#### HVGs ####
message(Sys.time(), " - Running modelGeneVar()")
# fit mean-variance relationship
dec <- modelGeneVar(spe)

# visualize mean-variance relationship
fit <- metadata(dec)

pdf(here(plot_dir, "spe_mean_var_curve.pdf"))
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
dev.off()

# select top HVGs
top.hvgs <- c(map(c(p1 = 0.1, p2 = 0.2, p5 = 0.5), ~getTopHVGs(dec, prop = .x)),
              map(c(fdr5 = 0.05, fdr1 = 0.01), ~getTopHVGs(dec, fdr.threshold = .x)))

(top.hvgs.length <- map_int(top.hvgs, length))
save(top.hvgs, file = here(data_dir, "top_hvgs.Rdata"))

#### Reduced Dims ####
num_reduced_dims <- 50

# compute PCA
message(Sys.time(), " - Compute PCA")

## Run PCA with selected set of HVGs
walk2(top.hvgs[1:2], names(top.hvgs[1:2]), function(hvgs, hvg_name){
    
    message(Sys.time(), " - PCA ", hvg_name)
    spe <<- runPCA(spe, 
                   subset_row = hvgs,
                   ncomponents = num_reduced_dims,
                   name = paste0("PCA_", hvg_name))
    
})

pca_sample <- plotReducedDim(spe, dimred = "PCA_p1",
                             color_by = "sample_id") +
    theme_bw()
ggsave(pca_sample, filename = here(plot_dir, "PCA_sample.png"), width = 10)

precent_var <- map_dfr(names(top.hvgs[1:2]), ~data.frame(hvg = .x,
                                                         PC = seq(ncol(reducedDim(spe, paste0("PCA_", .x)))),
                                                         precentVar = attr(reducedDim(spe, paste0("PCA_", .x)), "percentVar")))

## plot PCA elbow
pca_elbow <- ggplot(precent_var, aes(PC, precentVar, color = hvg)) +
    geom_line() +
    theme_bw()

ggsave(pca_elbow, filename = here(plot_dir, "PCA_elbow.png"))

percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
message("PCA elbow: ", chosen.elbow)

## TSNE & UMAP
message(Sys.time(), " - TSNE")
spe <-
    runTSNE(spe,
            dimred = "PCA_p1",
            name = "TSNE"
    )
Sys.time()

message(Sys.time(), " - UMAP")
spe <- runUMAP(spe, 
               dimred = "PCA_p1")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()

####   Compute GLM-PCA ####
message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial", sorted = FALSE, batch = as.factor(spe$sample_id))
spe <- devianceFeatureSelection(spe, assay = "counts", fam = "poisson", sorted = FALSE, batch = as.factor(spe$sample_id)) 

## plot binomial and poison deviance in first 100 selected genes ?

## calculate residuals from binomial model

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals( # default params
    spe,
    assay = "counts",
    fam = "binomial",
    type = "deviance"
    # batch = as.factor(spe$sample_id)
)

## Get HVDG
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
        ncomponents = num_red_dims,
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

# slurmjobs::job_single('01_preprocess_Harmony', create_shell = TRUE, memory = '25G', command = "Rscript 01_preprocess_Harmony.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()