
library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scran")
library("scater")
library("harmony")

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
    # BPPARAM = MulticoreParam(4),
    block = spe$sample_id #,
    # block.BPPARAM = MulticoreParam(4)
)

table(spe$scran_quick_cluster)

## plot quick clusters?

## compute size factors
message(Sys.time(), " - Run computeSumFactors()")

spe <- computeSumFactors(spe,
                         clusters = spe$scran_quick_cluster ,
                         BPPARAM = MulticoreParam(4)
)
Sys.time()


message("summary sizeFactors")
summary(sizeFactors(spe))

spe <- logNormCounts(spe)

#### HVGs ####

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
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)

#### Reduced Dims ####

# compute PCA
message(Sys.time(), " - Compute PCA")
set.seed(20240618)
spe <- runPCA(spe, subset_row = top_hvgs)

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
message("PCA elbow: ", chosen.elbow)

pdf(file.path(plot_dir, "pca_elbow.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

## TSNE & UMAP
message(Sys.time(), " - TSNE")
set.seed(20240618)
spe <-
    runTSNE(spe,
            dimred = "PCA",
            name = "TSNE"
    )
Sys.time()

message(Sys.time(), " - UMAP")
set.seed(20240618)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()


#### Harmony Batch Correction ####
message(Sys.time(), " - Harmony")
set.seed(20240618)
spe <- RunHarmony(spe, "sample_id")
Sys.time()


message(Sys.time(), " - Harmony UMAP")
set.seed(20240618)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
Sys.time()
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


message(Sys.time(), " - HARMONY TSNE")
set.seed(20240618)
spe <-
    runTSNE(spe,
            dimred = "HARMONY",
            name = "TSNE.HARMONY"
    )
Sys.time()

saveRDS(spe, file = here("processed-data", "02_build_spe", "spe_harmony.rds"))

# slurmjobs::job_single('01_preprocess_Harmony', create_shell = TRUE, memory = '25G', command = "Rscript 01_preprocess_Harmony.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()