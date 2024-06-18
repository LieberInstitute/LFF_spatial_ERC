
library("ggplot2")
library("SpotSweeper")
library("spatialLIBD")
library("here")
library("sessioninfo")
library("ggpubr")


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

# ## Add Scuttle QC metrics
# # change from gene id to gene names
# rownames(spe) <- rowData(spe)$gene_name
# 
# # identifying the mitochondrial transcripts
# is.mito <- rownames(spe)[grepl("^MT-", rownames(spe))]
# 
# # calculating QC metrics for each spot using scuttle
# spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(Mito = is.mito))
# colData(spe)
# 
# ## subsets_Mito_percent == expr_chrM_ratio*100
# head(colData(spe)[,c("expr_chrM_ratio", "subsets_Mito_percent")], 50)
# mito_diff = spe$expr_chrM_ratio*100 - spe$subsets_Mito_percent
# summary(mito_diff)
# head(mito_diff)

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


table(spe$sum_umi_outliers, spe$scran_low_lib_size)
#         TRUE  FALSE
# FALSE   5994 116432
# TRUE     223     97
table(spe$sum_umi_outliers, spe$scran_low_lib_size_edge)
table(spe$sum_umi_outliers, spe$scran_low_lib_size_edge)

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

table(spe$scran_discard)
table(spe$local_outliers)
# FALSE   TRUE 
# 122343    403
table(spe$local_outliers, spe$edge_spot)


#### find artifacts using SpotSweeper ####
## Only works one sample at a time
# spe <- findArtifacts(spe,
#                      mito_percent = "expr_chrM_ratio",
#                      mito_sum = "expr_chrM",
#                      n_rings = 5,
#                      name = "artifact"
# )

# Error in .set_internal_all(x, value, getfun = int_colData, setfun = `int_colData<-`,  : 
#                                invalid 'value' in 'reducedDims(<SpatialExperiment>) <- value'
#                            each element of 'value' should have number of rows equal to 'ncol(x)'

set.seed(20240617)

artifact_df <- purrr::map_dfr(unique(spe$sample_id), function(samp){
    message(Sys.time(), " - ", samp)
    spe_temp <- findArtifacts(spe[,spe$sample_id == samp],
                              mito_percent = "expr_chrM_ratio",
                              mito_sum = "expr_chrM",
                              n_rings = 5,
                              name = "artifact"
    )
    
    # print(plotQC(spe_temp,
    #        metric = "expr_chrM_ratio",
    #        outliers = "artifact", 
    #        point_size = 1.1, 
    #        stroke = 0.75
    # ))
    
    return(as.data.frame(colData(spe_temp)[,c("sample_id", "key", "artifact")]))
})

head(artifact_df)

# colData(spe) <- merge(colData(spe), artifact_df)
identical(spe$key, artifact_df$key)
spe$artifact <- artifact_df$artifact


#### save spot sweeper data ####
spotsweeper_data <- as.data.frame(colData(spe)[,c("sample_id", "key", "array_row", "array_col", "sum_umi_outliers", "sum_gene_outliers", "expr_chrM_ratio_outliers", "local_outliers", "artifact")])
head(spotsweeper_data)
write.csv(spotsweeper_data, file = here(data_dir, "SpotSweeper_data.csv"))

# dryspots <- c("V13B23-364_A1", "V13B23-365_B1", "V13B23-363_B1")
# Br5599 - artifact is not in corner (not sure this is a dry corner)
# Br6085 - finds vascular tissue 
# Br6098 - flagged corner + top edge, but lots of other tissue
table(spe$BrNum[which(spe$Visium_slide %in% dryspots)])

#### plotting ####
plot_all_spot_sweep <- function(spe, sample = unique(spe$sample_id)[1]){
    
    spe <- spe[,spe$sample_id == sample]
    
    # library size
    p1 <- plotQC(spe,
                 metric = "sum_umi_log",
                 outliers = "sum_umi_outliers", point_size = 1.1
    ) +
        ggtitle(paste(sample, "Sum UMI"))
    
    # unique genes
    p2 <- plotQC(spe,
                 metric = "sum_gene_log",
                 outliers = "sum_gene_outliers", point_size = 1.1
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
                 metric = "sum_umi_log",
                 outliers = "local_outliers", point_size = 1.1, stroke = 0.75
    ) +
        ggtitle("All Local Outliers")
    
    ## artifact
    p5 <- plotQC(spe,
                 metric = "sum_umi_log",
                 outliers = "artifact", point_size = 1.1, stroke = 0.75
    ) +
        ggtitle("Artifact")
    
    # plot
    plot_list <- list(p1, p2, p3, p4, p5)
    ggarrange(
        plotlist = plot_list,
        ncol = 3, nrow = 2,
        common.legend = FALSE
    )
} 

pdf(here(plot_dir, "SpotSweeper_ERC.pdf"), width = 12, height = 10)
purrr::map(sort(unique(spe$sample_id)), ~plot_all_spot_sweep(spe = spe, sample = .x))
dev.off()

## Consistent Outliers
spotsweeper_data |>
    as.tibble() |>
    filter(local_outliers) |>
    # count(array_row, array_col, sum_umi_outliers, sum_gene_outliers, expr_chrM_ratio_outliers) |>
    count(array_row, array_col) |>
    arrange(-n)
    
#   array_row array_col     n
# 1        47        77    28
# 2        38        88    26
# 3        57        99    24
# 4        19        91    18
# 5        49       103    15
# 6        28        32    11
# 7        53       117    10


