
library("spatialLIBD")
library("tidyverse")
library("scran") 
library("scater")
library("here")
library("sessioninfo")
library("patchwork")
library("ggridges")

#### set up dirs, load spe ####
plot_dir <- here("plots", "02_build_spe", "05_QC_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe", "05_QC_spe")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load spe
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))

## fix scran_high_subsets_Mito_percent colname
colnames(colData(spe))[grep("scran_high_subsets_Mito_percent", colnames(colData(spe)))] <- "scran_high_Mito_percent"

#### In Tissue ####
message(Sys.time(), "- Explore Metrics In and Out Tissue")

sample_in_tissue <- table(spe$sample_id[spe$in_tissue])
sort(sample_in_tissue)
# Br2582 Br3974 Br5367 Br2305 Br6263 Br5712 Br1556 Br5832 Br5854 Br6321 Br1289 Br1706 Br5517 Br5415 Br1039 Br6476 Br6085 Br5276 Br5460 
# 2396   2728   2843   3082   3231   3406   3534   3560   3580   3647   3663   3754   3851   3865   3919   3989   4031   4057   4059 
# Br1691 Br6538 Br6423 Br5529 Br6098 Br5212 Br5161 Br5599 Br5426 Br5634 Br5941 Br6161 
# 4186   4214   4335   4583   4601   4677   4700   4820   4840   4843   4872   4880

summary(as.integer(sample_in_tissue))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2396    3570    3989    3960    4592    4880

sample_order = sort(unique(spe$sample_id))

vis_grid_clus(
    spe = spe,
    clustervar = "overlaps_tissue",
    pdf = here::here(plot_dir, "spe_erc_grid-in_tissue.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c("in" = "grey75", "out" = "orange"),
    sample_order = sample_order
)


#### UMI ####
vis_grid_gene(
    spe = spe,
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_umi_all.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

summary(spe$sum_umi[spe$in_tissue])

vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_umi_in_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

summary(spe$sum_umi[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1770    3327    3938    5334   39355 


vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_umi_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

summary(spe$sum_umi[!spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0   179.0   277.0   352.4   417.0  5625.0 


#### Mito Rate ####
vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "expr_chrM_ratio",
    pdf = here::here(plot_dir, "spe_erc_grid-expr_chrM_ratio_in_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

summary(spe$expr_chrM_ratio[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1747  0.2284  0.2587  0.3057  1.0000 


vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "expr_chrM_ratio",
    pdf = here::here(plot_dir, "spe_erc_grid-expr_chrM_ratio_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)


#### Sum Gene ####
vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "sum_gene",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_gene_in_tissue.pdf"),
    assayname = "counts",
    point_size = 1
)

summary(spe$sum_gene[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1     876    1573    1708    2361    7735

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "sum_gene",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_gene_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1
)

summary(spe$sum_gene[!spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0   121.0   182.0   231.6   277.0  2687.0

#### marker genes ####

vis_gene_MBP_example <- vis_gene(
    spe = spe[,spe$in_tissue],
    sampleid = "Br5460",
    geneid = "MBP",
    assayname = "counts",
    point_size = 1.7
)

ggsave(vis_gene_MBP_example, filename = here(plot_dir, "vis_gene_MBP_Br5460.png"))

vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "MBP",
    pdf = here::here(plot_dir, "spe_erc_grid_gene-MBP.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

vis_grid_clus(
    spe = spe[,spe$in_tissue],
    clustervar = "10x_kmeans_9_clusters",
    pdf = here::here(plot_dir, "spe_erc_grid_clus-10x_kmeans_9_clusters.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

#### Outliers ####
message(Sys.time(), "- Scran Outlier Violin Plots")

pdf(here(plot_dir, "spe_erc_QC_outliers_violin_plots.pdf"), width = 21)

plotColData(spe[,spe$in_tissue], x = "sample_id", y = "sum_umi", colour_by = "scran_low_lib_size") +
    ggtitle("sum_umi") +
    facet_wrap(~ spe$round[spe$in_tissue], scales = "free_x", nrow = 2)

plotColData(spe[,spe$in_tissue], x = "sample_id", y = "sum_gene", colour_by = "scran_low_n_features") +
    ggtitle("sum_gene") +
    facet_wrap(~ spe$round[spe$in_tissue], scales = "free_x", nrow = 2)

plotColData(spe[,spe$in_tissue], x = "sample_id", y = "expr_chrM_ratio", colour_by = "scran_high_Mito_percent") +
    ggtitle("expr_chrM_ratio") +
    facet_wrap(~ spe$round[spe$in_tissue], scales = "free_x", nrow = 2)

dev.off()

vis_grid_clus(
    spe = spe,
    clustervar = "scran_low_lib_size",
    pdf = here::here(plot_dir, "spe_erc_grid_qc-low_lib_size.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)

vis_grid_clus(
    spe = spe,
    clustervar = "scran_low_n_features",
    pdf = here::here(plot_dir, "spe_erc_grid_qc-low_n_features.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)

vis_grid_clus(
    spe = spe,
    clustervar = "scran_high_Mito_percent",
    pdf = here::here(plot_dir, "spe_erc_grid_qc-high_Mito_percent.pdf"),
    sort_clust = FALSE,
    point_size = 1,
    colors = c(`TRUE` = "orange", `FALSE` = "lightskyblue3")
)

pd <- as.data.frame(colData(spe))

message("Scran UMI cutoffs:")
pd |> 
    filter(in_tissue, as.logical(scran_low_lib_size)) |> 
    group_by(sample_id) |> 
    dplyr::summarize(umi_cutoff = max(sum_umi)) |>
    arrange(umi_cutoff)

# sample_id umi_cutoff
# <chr>          <dbl>
# 1 Br6538            53
# 2 Br3974            75
# 3 Br5276            82
# 4 Br6321           224
# 5 Br1039           267
# 6 Br5529           294
# 7 Br6263           416
# 8 Br2582           441
# 9 Br5832           459
# 10 Br5517           482

#### ggplots for qc ####
median_qc_metrics <- pd |> 
    filter(in_tissue) |>
    select(sample_id, round, sum_umi, sum_gene, expr_chrM_ratio) |>
    pivot_longer(!c(sample_id, round), names_to = "metric", values_to = "value") |>
    group_by(sample_id, round, metric) |>
    summarise(median = median(value))

median_qc_metrics_boxplot <- median_qc_metrics |>
    ggplot(aes(x = round, y = median)) +
    geom_boxplot(outlier.shape = NULL) +
    geom_point() +
    facet_wrap(~metric, scales = "free_y")
    
ggsave(median_qc_metrics_boxplot, filename = here(plot_dir, "median_qc_metrics_boxplot.png"), height = 5)

# pd |>
#     filter(in_tissue) |>
#     count(sample_id) |>
#     ggplot(aes(y = n)) +
#     geom_histogram()

in_tissue_bar_preQC <- pd |>
    filter(in_tissue) |>
    count(sample_id) |>
    ggplot(aes(x= reorder(sample_id, n), y = n)) +
    geom_col() +
    coord_flip() +
    labs(title = "Pre-QC", x = "Sample", y= "n in-tissue spots") +
    theme_bw()
    
ggsave(in_tissue_bar_preQC, filename = here(plot_dir, "in_tissue_bar_preQC.png"), width =5)

in_tissue_boxplot_preQC <- pd |>
    filter(in_tissue) |>
    count(APOE, sample_id) |>
    ggplot(aes(x= APOE, y = n)) +
    geom_boxplot() +
    geom_point() +
    geom_text_repel(aes(label = sample_id), color = "grey50") +
    labs(title = "Pre-QC", y= "n in-tissue spots") +
    theme_bw()

ggsave(in_tissue_boxplot_preQC, filename = here(plot_dir, "in_tissue_boxplot_preQC.png"), width =7)


message(Sys.time(), "- ggridge plots")


#### ggridge plots ####
qc_ggridge_sum_umi <- pd |>
    ggplot(aes(x = sum_umi, y= sample_id, fill = in_tissue, color = in_tissue)) +
    geom_density_ridges(quantile_lines = TRUE,
                        alpha = .5)  +
    scale_fill_manual(values = c(`TRUE` = "grey20", `FALSE` = "darkorange")) +
    scale_color_manual(values = c(`TRUE` = "grey20", `FALSE` = "darkorange"))+
    theme(legend.position = "bottom")

ggsave(qc_ggridge_sum_umi, filename = here(plot_dir, "spe_erc_ggridge-sum_umi.png"), height = 12)

qc_ggridge_sum_gene <- pd |>
    ggplot(aes(x = sum_gene, y= sample_id, fill = in_tissue, color = in_tissue)) +
    geom_density_ridges(quantile_lines = TRUE,
                        alpha = .5)  +
    scale_fill_manual(values = c(`TRUE` = "grey20", `FALSE` = "darkorange")) +
    scale_color_manual(values = c(`TRUE` = "grey20", `FALSE` = "darkorange"))+
    theme(legend.position = "bottom")

ggsave(qc_ggridge_sum_gene, filename = here(plot_dir, "spe_erc_ggridge-sum_gene.png"), height = 12)

qc_ggridge_expr_chrM_ratio <- pd |>
    ggplot(aes(x = expr_chrM_ratio, y= sample_id, fill = in_tissue, color = in_tissue)) +
    geom_density_ridges(quantile_lines = TRUE,
                        alpha = .5)  +
    scale_fill_manual(values = c(`TRUE` = "grey20", `FALSE` = "darkorange")) +
    scale_color_manual(values = c(`TRUE` = "grey20", `FALSE` = "darkorange")) +
    theme(legend.position = "bottom")

ggsave(qc_ggridge_expr_chrM_ratio, filename = here(plot_dir, "spe_erc_ggridge-expr_chrM_ratio.png"), height = 12)


ggsave(qc_ggridge_sum_umi + theme(legend.position = "None")+ 
           qc_ggridge_sum_gene + theme(axis.title.y=element_blank(),axis.text.y=element_blank(), legend.position = "bottom") +
           qc_ggridge_expr_chrM_ratio  + theme(axis.title.y=element_blank(),axis.text.y=element_blank(), legend.position = "None"), 
       filename = here(plot_dir, "spe_erc_ggridge.png"), 
       height = 12, width = 12)

## in tissue
qc_ggridge_sum_umi_in <- pd |>
    filter(in_tissue) |>
    ggplot(aes(x = sum_umi, y= sample_id)) +
    geom_density_ridges(quantile_lines = TRUE,
                        alpha = .5)  

ggsave(qc_ggridge_sum_umi_in, filename = here(plot_dir, "spe_erc_ggridge_in-sum_umi.png"), height = 12)

qc_ggridge_sum_gene_in <- pd |>
    filter(in_tissue) |>
    ggplot(aes(x = sum_gene, y= sample_id)) +
    geom_density_ridges(quantile_lines = TRUE,
                        alpha = .5)  

ggsave(qc_ggridge_sum_gene_in, filename = here(plot_dir, "spe_erc_ggridge_in-sum_gene.png"), height = 12)

qc_ggridge_expr_chrM_ratio_in <- pd |>
    filter(in_tissue) |>
    ggplot(aes(x = expr_chrM_ratio, y= sample_id)) +
    geom_density_ridges(quantile_lines = TRUE,
                        alpha = .5)  

ggsave(qc_ggridge_expr_chrM_ratio_in, filename = here(plot_dir, "spe_erc_ggridge_in-expr_chrM_ratio_umi.png"), height = 12)

ggsave(qc_ggridge_sum_umi_in +
           qc_ggridge_sum_gene_in + 
           qc_ggridge_expr_chrM_ratio_in,
       filename = here(plot_dir, "spe_erc_ggridge_in.png"), 
       height = 12, width = 12)


#### Spot Sweeper Data ####
message(Sys.time(), "- Add Spot Sweeper Data")

spotsweeper_data <- read.csv(here("processed-data", "02_build_spe", "04_SpotSweeper", "SpotSweeper_data.csv"), row.names = 1)

pd <- pd |>
    left_join(spotsweeper_data) 


pd |> 
    filter(in_tissue) |> 
    count(scran_discard, local_outliers) |>
    mutate(precent = 100*n/sum(n))
# scran_discard local_outliers      n     precent
# 1          TRUE          FALSE   6663  5.42828280
# 2          TRUE           TRUE    285  0.23218679
# 3         FALSE          FALSE 115680 94.24339693
# 4         FALSE           TRUE    118  0.09613348

#### Read in Annotations ####
message(Sys.time(), "- Add Manual annotations")

qc_anno <- map_dfr(list.files(here("processed-data", "02_build_spe", "02_QC_app"), full.names = TRUE), 
                   ~read.csv(.x) |> 
                       mutate(file = gsub("spatialLIBD_ManualAnnotation_|.csv", "",basename(.x)))) |>
    filter(!is.na(ManualAnnotation),
           ManualAnnotation != "qc_drop", ## qc_drop was an old annotation
           spot_name != "TAACATACAATGTGGG-1_Br5367" ## mis-classified - fine
    ) |> 
    mutate(qc_anno = case_when(
        grepl("fold", ManualAnnotation) ~ "fold",
        grepl("out_edge", ManualAnnotation) ~ "out_edge",
        grepl("tissue", ManualAnnotation) ~ "out_tissue",
        grepl("torn_edge|circle", ManualAnnotation) ~ "vessel",
        TRUE ~ ManualAnnotation)
    )

qc_anno |> dplyr::count(file)
qc_anno |> dplyr::count(ManualAnnotation, qc_anno)
# qc_anno |> 
#     filter(ManualAnnotation == "qc_drop") |>
#     group_by(sample_id, file) |>
#     summarize(n = n(),
#               anno = paste(unique(ManualAnnotation), collapse = ",")
#               # ,
#               # files = paste(unique(file), collapse = ",")
#               )

qc_anno_clean <- qc_anno |>
    group_by(sample_id, spot_name) |>
    arrange(desc(file)) |>
    slice(1) |> ## select most recent anotations if there are more than one for a spot
    select(sample_id, sample_id, spot_name, qc_anno, file)  |>
    left_join(pd |>
                  select(sample_id, spot_name = key, in_tissue, scran_low_lib_size_edge, scran_discard)) |>
    filter(in_tissue) |>
    as.data.frame()

qc_anno_clean |>
    count(qc_anno) |>
    mutate(precent = n/sum(spe$in_tissue))
#      qc_anno    n
# 1       None    2
# 2       fold 1419
# 3   out_edge  257
# 4 out_tissue   69
# 5       tail   74
# 6     vessel  944

qc_anno_clean |>
    count(sample_id, qc_anno, scran_discard)

qc_anno_clean |>
    count(qc_anno, scran_discard)

## spots with more than 1 annotation
stopifnot(all(qc_anno_clean  |> group_by(spot_name) |> count() |> pull(n) == 1))

## test
# qc_anno_clean  |> filter(spot_name == "AACTAGCGTATCGCAC-1_Br5517")
# qc_anno_clean  |> filter(sample_id == "Br5276") |> arrange(file)


write.csv(qc_anno_clean, file = here(data_dir, "spe_qc_anno_clean.csv"), row.names = FALSE)

# qc_anno |> filter((spot_name == "TCTTTCCTTCGAGATA-1_Br5367" & ManualAnnotation == "out_tissue"))

## Add qc_anno to pd
# pd$qc_anno <- NULL ## reset
pd <- pd |>
    left_join(qc_anno_clean |> select(qc_anno, key = spot_name)) |>
    mutate(qc_anno = factor(ifelse(is.na(qc_anno) & in_tissue, "None", qc_anno)),
           scran_qc_anno = factor(ifelse(as.logical(scran_discard), paste0("scran-", qc_anno), as.character(qc_anno))),
           ss_qc_anno = factor(ifelse(local_outliers, paste0("ss-", qc_anno), as.character(qc_anno))),
           qc_anno_all = case_when(local_outliers ~ ss_qc_anno,
                                   as.logical(scran_low_lib_size_edge) & qc_anno == "None" ~ "scran-edge",
                                   as.logical(scran_discard) ~ scran_qc_anno,
                                   TRUE ~ qc_anno
                                   ),
           drop = local_outliers | (qc_anno %in% c("out_tissue", "tail"))
    )

# qc_anno_clean|> count(qc_anno)
pd|> count(in_tissue, qc_anno)
pd|> count(in_tissue, qc_anno, scran_qc_anno)
pd |> count(in_tissue, qc_anno_all)

pd |> 
    count(scran_discard) |>
    mutate(percent = 100*n/ncol(spe))

pd |> 
    filter(in_tissue) |>
    count(drop) |>
    mutate(percent = 100*n/sum(spe$in_tissue))

# drop      n    percent
# 1 FALSE 122202 99.5568084
# 2  TRUE    544  0.4431916

table(pd$qc_anno, pd$local_outliers)

## add to spe
spe$local_outliers <- pd$local_outliers
spe$qc_anno <- pd$qc_anno
spe$scran_qc_anno <- pd$scran_qc_anno
spe$ss_qc_anno <- pd$ss_qc_anno
spe$qc_anno_all <- pd$qc_anno_all

table(spe$qc_anno)
# fold       None   out_edge out_tissue       tail     vessel 
# 1419     119983        257         69         74        944 
table(spe$scran_qc_anno)
table(spe$qc_anno_all)


qc_colors <- c(`ss-out_edge` = "#06F8EC",
               `scran-out_edge` = "#87BBFF",
               out_edge = "#097FE0",
               `ss-fold` = "#91FA11",
               `scran-fold` = "#05F8A3",
               fold = "#4BA402",
               `ss-out_tissue` ="#FFCAD4",
               `scran-out_tissue` ="#FB9D88",
               out_tissue ="#F24018",
               `ss-tail` = "#FFF899",
               `scran-tail` = "#EDDD00",
               tail = "#FE7215",
               `ss-vessel` ="#FD6092",
               `scran-vessel` ="#FA73E6",
               vessel ="#BB1EF4",
               `scran-edge` = "brown",
               `ss-None` = "white",
               `scran-None` = "black",
               None = "grey50")

#### bar plots ####

pd |> 
    filter(in_tissue) |>
    ggplot(aes(x = sample_id, fill = qc_anno_all)) +
    geom_bar() +
    scale_fill_manual(values = qc_colors) +
    coord_flip()



pdf(here(plot_dir, "spe_erc_QC_annotations.pdf"))
map(sample_order, 
    ~vis_clus(
        spe = spe,
        sampleid = .x,
        point_size = 1.2,
        clustervar = "qc_anno_all",
        colors = qc_colors
    ))
dev.off()

# summary
# qc_summary <- qc_anno_clean |> 
#     count(sample_id, qc_anno) |>
#     arrange(qc_anno) |>
#     mutate(summary = paste0(qc_anno, ":" ,n)) |>
#     group_by(sample_id) |>
#     summarize(n = sum(n),
#               anno = paste(summary, collapse = ","))

qc_summary <- pd |> 
    filter(in_tissue) |>
    count(sample_id, qc_anno, local_outliers, drop) |>
    group_by(sample_id) |>
    mutate(precent = 100*n/sum(n))

write.csv(qc_summary, file = here(data_dir, "spe_qc_summary.csv"), row.names = FALSE)


pd |> 
    filter(in_tissue,) |>
    count(sample_id, drop) |>
    group_by(sample_id) |>
    mutate(precent = 100*n/sum(n)) |>
    filter(drop) |>
    arrange(-precent)

# sample_id drop      n precent
# <chr>     <lgl> <int>   <dbl>
#     1 Br1706    TRUE     55   1.47 
# 2 Br2582    TRUE     29   1.21 
# 3 Br5276    TRUE     46   1.13 
# 4 Br5712    TRUE     33   0.969
# 5 Br5517    TRUE     28   0.727
    
#### Drop spots ####
message(Sys.time(), "- Drop Spots")

pd |> count(drop)
# Drop      n
# 1 FALSE 122202
# 2  TRUE    544
# 3    NA  32006

pd |> count(qc_anno, local_outliers, drop)
# qc_anno local_outliers  Drop      n
# 1        fold          FALSE FALSE   1408
# 2        fold           TRUE  TRUE     11
# 3        None          FALSE FALSE 119633
# 4        None           TRUE  TRUE    350
# 5    out_edge          FALSE FALSE    221
# 6    out_edge           TRUE  TRUE     36
# 7  out_tissue          FALSE  TRUE     67
# 8  out_tissue           TRUE  TRUE      2
# 9        tail          FALSE  TRUE     74
# 10     vessel          FALSE FALSE    940
# 11     vessel           TRUE  TRUE      4
# 12       <NA>             NA FALSE  32006


pd_qc <- pd |> 
    filter(in_tissue, !drop) 

table(pd_qc$sample_id)
# Br1039 Br1289 Br1556 Br1691 Br1706 Br2305 Br2582 Br3974 Br5161 Br5212 Br5276 Br5367 Br5415 Br5426 Br5460 Br5517 Br5529 Br5599 Br5634 Br5712 
# 3906   3648   3524   4161   3699   3069   2367   2722   4682   4674   4011   2828   3858   4824   4032   3823   4570   4812   4830   3373 
# Br5832 Br5854 Br5941 Br6085 Br6098 Br6161 Br6263 Br6321 Br6423 Br6476 Br6538 
# 3546   3575   4869   4017   4570   4864   3226   3641   4315   3964   4202 

pd_qc |> count(sample_id) |> summary()

median_PostQC_metrics <- pd_qc |> 
    select(sample_id, round, sum_umi, sum_gene, expr_chrM_ratio) |>
    pivot_longer(!c(sample_id, round), names_to = "metric", values_to = "value") |>
    group_by(sample_id, round, metric) |>
    summarise(median_postQC = median(value)) |>
    left_join(median_qc_metrics |> rename(median_preQC = median))

median_PostQC_metrics |>
    mutate(diff = median_postQC - median_preQC) |>
    group_by(metric) |>
    summarise(mean(diff))

median_PostQC_metrics_boxplot <- median_PostQC_metrics |>
    pivot_longer(!c(sample_id, round, metric), names_to = "QC", values_to = "median", names_prefix = "median_") |>
    mutate(QC = factor(QC, levels = c("preQC", "postQC"))) |>
    ggplot(aes(x = QC, y = median)) +
    geom_boxplot(outlier.shape = NULL) +
    geom_point() +
    geom_line(aes(group = sample_id), color = "skyblue") +
    facet_wrap(~metric, scales = "free_y")

ggsave(median_PostQC_metrics_boxplot, filename = here(plot_dir, "median_PostQC_metrics_boxplot.png"), height = 7, width = 9)


spe <- spe[,pd_qc$key]

message("QCed total spots: ", ncol(spe))

message("QCed spot summary: ")
sample_in_tissue <- table(spe$sample_id[spe$in_tissue])
sort(sample_in_tissue)

summary(as.integer(sample_in_tissue))

message(Sys.time(), "- Save QCed SPE")
saveRDS(spe, file = here("processed-data", "02_build_spe", "spe.rds"))

# slurmjobs::job_single('05_QC_spe', create_shell = TRUE, memory = '25G', command = "Rscript 05_QC_spe.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
