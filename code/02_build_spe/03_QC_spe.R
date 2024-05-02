
library("spatialLIBD")
library("tidyverse")
library("scran") 
library("scater")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "02_build_spe", "03_QC_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load spe
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))

spe$trimmed <- ifelse(grepl("untrimmed", spe$base_path), "untrimmed", "trimmed")

#### Tissue Overlap ####

(sample_in_tissue <- table(spe$sample_id[spe$in_tissue]))
# V13B23-363_A1 V13B23-363_B1 V13B23-363_C1 V13B23-363_D1 V13B23-364_A1 V13B23-364_B1 V13B23-364_C1 V13B23-364_D1 
# 4872          4820          3560          4843          4031          4880          3580          3663 
# V13B23-365_A1 V13B23-365_B1 V13B23-365_C1 V13B23-365_D1 V13B23-366_B1 V13B23-366_C1 V13B23-366_D1 V13Y24-340_A1 
# 4335          4601          4700          4583          3647          4214          3231          3851 
# V13Y24-340_B1 V13Y24-340_C1 V13Y24-340_D1 V13Y24-342_A1 V13Y24-342_B1 V13Y24-342_C1 V13Y24-342_D1 V13Y24-343_A1 
# 3082          3406          2843          3989          4840          3754          4059          4677 
# V13Y24-343_B1 V13Y24-343_C1 V13Y24-343_D1 V13Y24-344_A1 V13Y24-344_B1 V13Y24-344_C1 V13Y24-344_D1 
# 3919          3534          4057          2396          3865          2728          4186 

summary(as.integer(sample_in_tissue))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2396    3570    3989    3960    4592    4880

vis_grid_clus(
    spe = spe,
    clustervar = "overlaps_tissue",
    pdf = here::here(plot_dir, "in_tissue_grid.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c("in" = "grey90", "out" = "orange")
)

summary(spe$sum_umi[!inTissue(spe_raw)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.0   101.0   162.0   231.3   287.8  1553.0



#### UMI ####
summary(spe$sum_umi[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1770    3327    3938    5334   39355 

summary(spe$sum_umi[!spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0   179.0   277.0   352.4   417.0  5625.0 

## spot plots 
vis_grid_gene(
    spe = spe,
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "sum_umi_all.pdf"),
    # assayname = "counts",
    point_size = 1
)

vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "sum_umi_in_tissue.pdf"),
    # assayname = "counts",
    point_size = 1
)

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "sum_umi_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1
)

## Violin plots 


#### Mito Rate ####
summary(spe$expr_chrM_ratio[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1747  0.2284  0.2587  0.3057  1.0000 

summary(spe$expr_chrM_ratio[!spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1866  0.2444  0.2691  0.3169  0.8623 

vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "expr_chrM_ratio",
    pdf = here::here(plot_dir, "expr_chrM_ratio_in_tissue.pdf"),
    # assayname = "counts",
    point_size = 1
)

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "expr_chrM_ratio",
    pdf = here::here(plot_dir, "expr_chrM_ratio_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1
)


#### Sum Gene ####
summary(spe$sum_gene[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1     876    1573    1708    2361    7735

summary(spe$sum_gene[!spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0   121.0   182.0   231.6   277.0  2687.0

vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "sum_gene",
    pdf = here::here(plot_dir, "sum_gene_in_tissue.pdf"),
    point_size = 1
)

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "sum_gene",
    pdf = here::here(plot_dir, "sum_gene_out_tissue.pdf"),
    point_size = 1
)


#### Outliers ####

spe <- spe[, spe$in_tissue]

spe$low_lib_size = isOutlier(
    spe$sum_umi,
    type = "lower",
    log = TRUE,
    batch = spe$sample_id
)

spe$low_n_features = isOutlier(
    spe$sum_gene,
    type = "lower",
    log = TRUE,
    batch = spe$sample_id
)

spe$high_Mito_percent = isOutlier(
    spe$expr_chrM_ratio,
    type = "higher",
    batch = spe$sample_id
)

pdf(here(plot_dir, "QC_outliers.pdf"), width = 21)

plotColData(spe, x = "sample_id", y = "sum_umi", colour_by = "low_lib_size") +
    ggtitle("sum_umi") +
    facet_wrap(~ spe$trimmed[spe$in_tissue], scales = "free_x", nrow = 2)

plotColData(spe, x = "sample_id", y = "sum_gene", colour_by = "low_n_features") +
    ggtitle("sum_gene") +
    facet_wrap(~ spe$trimmed[spe$in_tissue], scales = "free_x", nrow = 2)

plotColData(spe, x = "sample_id", y = "expr_chrM_ratio", colour_by = "high_Mito_percent") +
    ggtitle("expr_chrM_ratio") +
    facet_wrap(~ spe$trimmed[spe$in_tissue], scales = "free_x", nrow = 2)

dev.off()

vis_grid_clus(
    spe = spe,
    clustervar = "low_lib_size",
    pdf = here::here(plot_dir, "QC_grid_low_lib_size.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)

vis_grid_clus(
    spe = spe,
    clustervar = "low_n_features",
    pdf = here::here(plot_dir, "QC_grid_low_n_features.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)

vis_grid_clus(
    spe = spe,
    clustervar = "high_Mito_percent",
    pdf = here::here(plot_dir, "QC_grid_high_Mito_percent.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)



## Scran QC
qcstats <- perCellQCMetrics(spe, subsets = list(
    Mito = which(seqnames(spe) == "chrM")
))
qcfilter <- quickPerCellQC(qcstats, sub.fields = "subsets_Mito_percent")
colSums(as.matrix(qcfilter))



## Find edge spots
spot_coords <- colData(spe) |>
    as.data.frame() |>
    select(sample_id, array_row, array_col)

edge_spots <- spot_coords |>
    group_by(sample_id, array_row) |>
    mutate(edge_col = array_col == min(array_col) | array_col == max(array_col)) |>
    group_by(sample_id, array_col) |>
    mutate(edge_row = array_row == min(array_row) | array_row == max(array_row)) |>
    group_by(sample_id) |>
    mutate(edge = edge_row | edge_col)

edge_spots |> count(edge)

spe$edge_spot <- edge_spots$edge

vis_grid_clus(
    spe = spe,
    clustervar = "edge_spot",
    pdf = here::here(plot_dir, "QC_edge.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)

