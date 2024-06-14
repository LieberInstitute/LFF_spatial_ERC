
library("SpotSweeper")
library("spatialLIBD")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "02_build_spe", "04_SpotSweeper")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe", "04_SpotSweeper")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))
colnames(colData(spe))

# drop out-of-tissue spots
table(spe$in_tissue)
spe <- spe[, spe$in_tissue]
dim(spe)
# [1]  30494 122746

## Add Scuttle QC metrics
# change from gene id to gene names
rownames(spe) <- rowData(spe)$gene_name

# identifying the mitochondrial transcripts
is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]

# calculating QC metrics for each spot using scuttle
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
colData(spe)

## subsets_Mito_percent == expr_chrM_ratio*100
head(colData(spe)[,c("expr_chrM_ratio", "subsets_Mito_percent")], 50)
mito_diff = spe$expr_chrM_ratio*100 - spe$subsets_Mito_percent
summary(mito_diff)
head(mito_diff)

#### Spot sweeper ####
# library size
spe <- localOutliers(spe,
                     metric = "sum_umi",
                     direction = "lower",
                     log = TRUE
)

# unique genes
spe <- localOutliers(spe,
                     metric = "sum_gene",
                     direction = "lower",
                     log = TRUE
)

# mitochondrial percent
spe <- localOutliers(spe,
                     metric = "expr_chrM_ratio",
                     direction = "higher",
                     log = FALSE
)

table(spe$sum_umi_outliers)
table(spe$sum_umi_outliers, spe$scran_low_lib_size)
#         TRUE  FALSE
# FALSE   5994 116432
# TRUE     223     97

table(spe$sum_gene_outliers)
table(spe$sum_gene_outliers, spe$scran_low_lib_size)

table(spe$expr_chrM_ratio_outliers)
# FALSE   TRUE 
# 122703     43
table(spe$expr_chrM_ratio_outliers, spe$scran_high_subsets_Mito_percent)
#        TRUE  FALSE
# FALSE    458 122245
# TRUE      11     32

# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_umi_outliers) |
    as.logical(spe$sum_gene_outliers) |
    as.logical(spe$expr_chrM_ratio_outliers)

library(ggpubr)

# library size
p1 <- plotQC(spe,
             metric = "sum_umi_log",
             outliers = "sum_umi_outliers", point_size = 1.1
) +
    ggtitle("Sum UMI")

# unique genes
p2 <- plotQC(spe,
             metric = "sum_gene_log",
             outliers = "detected_outliers", point_size = 1.1
) +
    ggtitle("Sum Genes")

# mitochondrial percent
p3 <- plotQC(spe,
             metric = "expr_chrM_ratio",
             outliers = "expr_chrM_ratio_outliers", point_size = 1.1
) +
    ggtitle("ChrM Ratio")

# all local outliers
p4 <- plotQC(spe,
             metric = "sum_log",
             outliers = "local_outliers", point_size = 1.1, stroke = 0.75
) +
    ggtitle("All Local Outliers")

# plot
plot_list <- list(p1, p2, p3, p4)
ggarrange(
    plotlist = plot_list,
    ncol = 2, nrow = 2,
    common.legend = FALSE
)

# find artifacts using SpotSweeper
spe <- findArtifacts(spe,
                     mito_percent = "expr_chrM_ratio",
                     mito_sum = "expr_chrM",
                     n_rings = 5,
                     name = "artifact"
)



