## Louise Huuki-Myers, April 2025
## Run differential proportion analysis with Crumblr

library("SingleCellExperiment")
library("tidyverse")
library("crumblr")
library("variancePartition")
library("here")
library("sessioninfo")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "22_crumblr_sn")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "22_crumblr_sn")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### Load data ####
load(here("processed-data", "04_snRNA-seq", "15_cluster_update_sce","cell_type_proportions.Rdata"), verbose = TRUE)
# cell_type_proportions

pd <- colData(sce_pb) |> as.data.frame()

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- cell_type_proportions |>
    mutate(error = sqrt((prop * (1-prop))/sum(n)),
           anno = paste0(n,"/", sum(n))) |>
    select(sample_id, cell_type_anno, error, anno)

## create cell count matrix
cell_counts <- cell_type_proportions |> 
    ungroup() |>
    select(sample_id, cell_type_anno, n) |>
    pivot_wider(names_from = cell_type_anno, values_from = n)  |>
    column_to_rownames("sample_id")

## replace NA w/ 0
cell_counts[is.na(cell_counts)] <- 0

## create info matrix

erc_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv")) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, exp_round, seq_round)|>
    column_to_rownames("sample_id")

setequal(rownames(erc_info), rownames(cell_counts))

erc_info <- erc_info[rownames(cell_counts), ]
identical(rownames(erc_info), rownames(cell_counts))

#### Apply crumblr transformation ####
cobj <- crumblr(cell_counts)

# E:numeric matrix of CLR transformed counts
clr_prop_long <- cobj$E |>
    as.data.frame() |>
    rownames_to_column("cell_type_anno") |>
    pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR") |>
    left_join(erc_info |> rownames_to_column("sample_id")) |>
    left_join(error)


# Partition variance into components for sample
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | exp_round)+ (1 | seq_round)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

ggsave(fig.vp, filename = here(plot_dir, "crumblr_cell_type_vp.png"))

#### PCA ####
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
#   shape by Stimulated vs unstimulated
ggplot(df_pca, aes(PC1, PC2, color = exp_round, shape = APOE_carrier)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

#### Hierarchical clustering ####
library("dendextend")

sce_pb2 <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_pseudobulk_only-cell_type_anno.rds"))
dim(sce_pb2)

dist.clusCollapsed <- dist(t(logcounts(sce_pb2)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)
# plot(dend, main = "hierarchical cluster dend", horiz = TRUE)

# Use variancePartition workflow to analyze each cell type
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance

# fit <- dream(cobj, ~ APOE_carrier + Sex + Age + Anc_Afr + exp_round + seq_round , erc_info)
fit <- dream(cobj, ~ APOE_carrier + Age + Anc_Afr + exp_round , erc_info)
# fit <- dream(cobj, ~ APOE_carrier, erc_info)
fit <- eBayes(fit)

# Extract results for each cell type
topTable(fit, coef = "APOE_carrierE4+", number = Inf)

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

# Plot hierarchy and testing results
plotTreeTest(res)

# Plot hierarchy and regression coefficients
plotTreeTestBeta(res)

combined_fig <- fig.vp +
    theme(legend.position = "left") |
    plotTreeTestBeta(res) +
    theme(legend.position = "bottom", legend.box = "vertical") |
    plotForest(res, hide = FALSE) 

ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined.png"))

# slurmjobs::job_single('22_crumblr_sn', create_shell = TRUE, memory = '25G', command = "Rscript 19_sn_heatmaps.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
