
library("ggplot2")
library("SpotSweeper")
library("spatialLIBD")
library("purrr")
library("here")
library("sessioninfo")
library("ggpubr")


plot_dir <- here("plots", "02_build_spe", "04_SpotSweeper")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe", "04_SpotSweeper")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


set.seed(20240617)
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))
colnames(colData(spe))

# drop out-of-tissue spots
table(spe$in_tissue)
spe <- spe[, spe$in_tissue]
dim(spe)
# [1]  30494 122746

#### Spot sweeper ####
# library size
message(Sys.time(), "Local Outliers - sum_umi")
spe <- localOutliers(spe,
                     metric = "sum_umi",
                     direction = "lower",
                     log = TRUE
)
table(spe$sum_umi_outliers)


# unique genes
message(Sys.time(), "Local Outliers - sum_gene")
spe <- localOutliers(spe,
                     metric = "sum_gene",
                     direction = "lower",
                     log = TRUE
)
table(spe$sum_gene_outliers)

# mitochondrial percent
message(Sys.time(), "Local Outliers - expr_chrM_ratio")
spe <- localOutliers(spe,
                     metric = "expr_chrM_ratio",
                     direction = "higher",
                     log = FALSE
)
table(spe$expr_chrM_ratio_outliers)

# combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_umi_outliers) |
    as.logical(spe$sum_gene_outliers) |
    as.logical(spe$expr_chrM_ratio_outliers)

message("Local Outliers")
table(spe$local_outliers)

message("Local Outliers on edge")
table(spe$local_outliers, spe$edge_spot)


#### find artifacts using SpotSweeper ####
## Only works one sample at a time
message(Sys.time(), " - findArtifact")
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

## Add artifact to spe
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
# table(spe$BrNum[which(spe$Visium_slide %in% dryspots)])

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

message(Sys.time(), " - Plot SpotSweeper QC for all samples")
pdf(here(plot_dir, "SpotSweeper_ERC.pdf"), width = 12, height = 10)
purrr::map(sort(unique(spe$sample_id)), ~plot_all_spot_sweep(spe = spe, sample = .x))
dev.off()

## Consistent Outliers
# spotsweeper_data |>
#     as.tibble() |>
#     filter(local_outliers) |>
#     # count(array_row, array_col, sum_umi_outliers, sum_gene_outliers, expr_chrM_ratio_outliers) |>
#     count(array_row, array_col) |>
#     arrange(-n)
    
#   array_row array_col     n
# 1        47        77    28
# 2        38        88    26
# 3        57        99    24
# 4        19        91    18
# 5        49       103    15
# 6        28        32    11
# 7        53       117    10

# slurmjobs::job_single('04_SpotSweeper', create_shell = TRUE, memory = '25G', command = "Rscript 04_SpotSweeper.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
