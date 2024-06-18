
library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scran")
library("scater")

plot_dir <- here("plots", "05_spe_correct_cluster", "01_preprocess_Harmony")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "02_build_spe", "01_preprocess_Harmony")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

spe <- readRDS(here("processed-data", "02_build_spe", "spe.rds"))

#### Logcounts ####
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000506 0.452028 0.845337 1.000000 1.352287 9.960555

# calculate logcounts and store in object
spe <- logNormCounts(spe)

#### HVGs ####

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

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
set.seed(20240618)

# compute PCA
spe <- runPCA(spe, subset_row = top_hvgs)
# compute UMAP
spe <- runUMAP(spe, dimred = "PCA")

reducedDimNames(spe)

ggplot(
    data.frame(
        reducedDim(
            spe, "PCA"
        )
    ),
    aes(x = PC1, y = PC2, color = spe$sum_umi)
) +
    geom_point(size = 1) +
    theme_bw()


ggplot(
    data.frame(
        reducedDim(
            spe, "UMAP"
        )
    ),
    aes(x = UMAP1, y = UMAP2, color = spe$sample_id)
) +
    geom_point(size = 1) +
    theme_bw()

umap <- data.frame(
    reducedDim(
        spe, "UMAP"
    )
) 
umap$sample_id <- spe$sample_id

umap |>
    filter(UMAP1 > 6) |>
    count(sample_id)




