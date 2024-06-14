
library("SpotSweeper")
library("spatialLIBD")
library("here")
library("sessioninfo")
library("escheR")

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

# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_outliers) |
    as.logical(spe$detected_outliers) |
    as.logical(spe$subsets_Mito_percent_outliers)


# find artifacts using SpotSweeper
spe <- findArtifacts(spe,
                     mito_percent = "expr_chrM_ratio",
                     mito_sum = "expr_chrM",
                     n_rings = 5,
                     name = "artifact"
)

