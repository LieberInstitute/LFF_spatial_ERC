## Louise Huuki-Myers, April 2025
## Run differential proportion analysis with Crumblr

library("SingleCellExperiment")
library("tidyverse")
library("crumblr")
library("variancePartition")
library("here")
library("sessioninfo")

## Prep directories
plot_dir <- here("plots", "05_spe_correct_cluster", "26_crumblr_SpD")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "05_spe_correct_cluster", "26_crumblr_SpD")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### Load data ####

sce_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds"))
dim(sce_pb)

pd <- colData(sce_pb) |> as.data.frame()

## create cell count matrix
cell_counts <- pd |> 
    select(sample_id, SpD_syn, ncells) |>
    pivot_wider(names_from = "SpD_syn", values_from = "ncells") |>
    column_to_rownames("sample_id")

## replace NA w/ 0
cell_counts[is.na(cell_counts)] <- 0

## create info matrix

erc_info <- pd |> select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, Visium_slide, round) |> unique()
rownames(erc_info) <- erc_info$sample_id

identical(rownames(erc_info), rownames(cell_counts))

#### Apply crumblr transformation ####
cobj <- crumblr(cell_counts)

# Partition variance into components for sample
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | Visium_slide)+ (1 | round)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

ggsave(fig.vp, filename = here(plot_dir, "crumblr_SpD_vp.png"))

#### PCA ####
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
#   shape by Stimulated vs unstimulated
ggplot(df_pca, aes(PC1, PC2, color = Visium_slide, shape = APOE_carrier)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

#### Hierarchical clustering ####
library("dendextend")

sce_pb2 <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_sn_pseudobulk","sce_pseudobulk_only-SpD.rds"))
dim(sce_pb2)

dist.clusCollapsed <- dist(t(logcounts(sce_pb2)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)
# plot(dend, main = "hierarchical cluster dend", horiz = TRUE)

# Use variancePartition workflow to analyze each cell type
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance

# fit <- dream(cobj, ~ APOE_carrier + Sex + Age + Anc_Afr + exp_round + seq_round , erc_info)
fit <- dream(cobj, ~ APOE_carrier + Age + Visium_slide, erc_info)
# fit <- dream(cobj, ~ APOE_carrier, erc_info)
fit <- eBayes(fit)

# Test APOE carrier
topTable(fit, coef = "APOE_carrierE4+", number = Inf)

## visium slide
topTable(fit, coef = 4:7, number = Inf)

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

# Plot hierarchy and testing results
plotTreeTest(res)
plotTreeTestBeta(res) 

# Plot hierarchy and regression coefficients
plotTreeTestBeta(res)

combined_fig <- fig.vp +
    theme(legend.position = "left") |
    plotTreeTestBeta(res) +
    theme(legend.position = "bottom", legend.box = "vertical") |
    plotForest(res, hide = FALSE) 

ggsave(combined_fig, filename = here(plot_dir, "crumblr_SpD_combined.png"))

# slurmjobs::job_single('26_crumblr_SpD', create_shell = TRUE, memory = '25G', command = "Rscript 26_crumblr_SpD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
