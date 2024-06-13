
library("spatialLIBD")
library("tidyverse")
library("scran") 
library("scater")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "02_build_spe", "03_QC_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe", "03_QC_spe")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load spe
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))


#### In Tissue ####
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
summary(spe$sum_umi[spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1770    3327    3938    5334   39355 

table(spe$sum_umi[spe$in_tissue] < 10)
table(spe$sample_id[spe$in_tissue], spe$sum_umi[spe$in_tissue] < 10)

summary(spe$sum_umi[!spe$in_tissue])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0   179.0   277.0   352.4   417.0  5625.0 

## spot plots 
vis_grid_gene(
    spe = spe,
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_umi_all.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

vis_grid_gene(
    spe = spe[,spe$in_tissue],
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_umi_in_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "sum_umi",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_umi_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

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
    pdf = here::here(plot_dir, "spe_erc_grid-expr_chrM_ratio_in_tissue.pdf"),
    # assayname = "counts",
    point_size = 1,
    sample_order = sample_order
)

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "expr_chrM_ratio",
    pdf = here::here(plot_dir, "spe_erc_grid-expr_chrM_ratio_out_tissue.pdf"),
    assayname = "counts",
    point_size = 1,
    sample_order = sample_order
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
    pdf = here::here(plot_dir, "spe_erc_grid-sum_gene_in_tissue.pdf"),
    point_size = 1
)

vis_grid_gene(
    spe = spe[,!spe$in_tissue],
    geneid = "sum_gene",
    pdf = here::here(plot_dir, "spe_erc_grid-sum_gene_out_tissue.pdf"),
    point_size = 1
)


#### Outliers ####
pdf(here(plot_dir, "spe_erc_QC_outliers_violin_plots.pdf"), width = 21)

plotColData(spe, x = "sample_id", y = "sum_umi", colour_by = "scran_low_lib_size") +
    ggtitle("sum_umi") +
    facet_wrap(~ spe$trimmed[spe$in_tissue], scales = "free_x", nrow = 2)

plotColData(spe, x = "sample_id", y = "sum_gene", colour_by = "scran_low_n_features") +
    ggtitle("sum_gene") +
    facet_wrap(~ spe$trimmed[spe$in_tissue], scales = "free_x", nrow = 2)

plotColData(spe, x = "sample_id", y = "expr_chrM_ratio", colour_by = "scran_high_Mito_percent") +
    ggtitle("expr_chrM_ratio") +
    facet_wrap(~ spe$trimmed[spe$in_tissue], scales = "free_x", nrow = 2)

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

## temp 

pdf(here(plot_dir, "spe_erc_Br5276-sum_umi.pdf"))
vis_gene(
    spe = spe[,spe$in_tissue],
    sampleid = "Br5276",
    geneid = "sum_umi",
    assayname = "counts"
    # colors = c("in" = "grey90", "out" = "orange")
)
dev.off()

pd <- as.data.frame(colData(spe))

pd |> 
    filter(in_tissue, as.logical(scran_low_lib_size)) |> 
    group_by(sample_id) |> 
    dplyr::summarize(umi_cutoff = max(sum_umi)) |>
    arrange(umi_cutoff)

## ggplots for qc

# qc_density_sum_umi <- pd |>
#     ggplot(aes(x = sum_umi, color = in_tissue)) +
#     geom_density() +
#     # facet_grid(sample_id~.)
#     facet_wrap(~sample_id, ncol = 2)
# 
# ggsave(qc_density_sum_umi, filename = here(plot_dir, "spe_erc_density-sum_umi.png"), height = 12)

library(ggridges)

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

library(patchwork)

ggsave(qc_ggridge_sum_umi + theme(legend.position = "None")+ 
           qc_ggridge_sum_gene + theme(axis.title.y=element_blank(),axis.text.y=element_blank(), legend.position = "bottom") +
           qc_ggridge_expr_chrM_ratio  + theme(axis.title.y=element_blank(),axis.text.y=element_blank(), legend.position = "None"), 
       filename = here(plot_dir, "spe_erc_ggridge.png"), 
       height = 12, width = 12)


#### Read in Annotations ####
qc_anno <- map_dfr(list.files(here("processed-data", "02_build_spe", "02_QC_app"), full.names = TRUE), 
                   ~read.csv(.x) |> mutate(file = gsub("spatialLIBD_ManualAnnotation_|.csv", "",basename(.x)))) |>
    filter(!is.na(ManualAnnotation),
           ManualAnnotation != "qc_drop", ## qc_drop was an old annotation
           spot_name != "TAACATACAATGTGGG-1_Br5367" ## mis-classified - fine
    ) |> 
    mutate(qc_anno = case_when(
        grepl("fold", ManualAnnotation) ~ "fold",
        grepl("out_edge", ManualAnnotation) ~ "out_edge",
        grepl("tissue", ManualAnnotation) ~ "out_tissue",
        TRUE ~ ManualAnnotation)
    )

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
    select(sample_id, sample_id, spot_name, qc_anno)  |>
    unique() |>
    left_join(pd |>
                  select(sample_id, spot_name = key, in_tissue, scran_low_lib_size_edge, scran_discard)) |>
    filter(in_tissue) |>
    as.data.frame()

qc_anno_clean |>
    count(qc_anno)
#      qc_anno    n
# 1     circle   45
# 2       fold 1255
# 3   out_edge  268
# 4 out_tissue   71
# 5       tail  105
# 6  torn_edge  820

qc_anno_clean |>
    count(sample_id, qc_anno, scran_discard)

qc_anno_clean |>
    count(qc_anno, scran_discard)

## spots with more than 1 annotation
qc_anno_clean  |> group_by(spot_name) |> filter(n() > 1) |> arrange(spot_name)

write.csv(qc_anno_clean, file = here(data_dir, "spe_qc_anno_clean.csv"), row.names = FALSE)

qc_summary <- qc_anno_clean |> 
    count(sample_id, qc_anno) |>
    arrange(qc_anno) |>
    mutate(summary = paste0(qc_anno, ":" ,n)) |>
    group_by(sample_id) |>
    summarize(n = sum(n),
              anno = paste(summary, collapse = ","))

write.csv(qc_summary, file = here(data_dir, "spe_qc_summary.csv"), row.names = FALSE)

# qc_anno |> filter((spot_name == "TCTTTCCTTCGAGATA-1_Br5367" & ManualAnnotation == "out_tissue"))
# 
# rownames(qc_anno_clean) <- qc_anno_clean$spot_name
# spe$qc_anno <- qc_anno_clean[colnames(spe),]$qc_anno
# spe$qc_anno[spe$in_tissue & is.na(spe$qc_anno)] <- "None"
# spe$qc_anno <- factor(spe$qc_anno)

pd <- pd |>
    left_join(qc_anno_clean |> select(qc_anno, key = spot_name)) |>
    mutate(qc_anno = factor(ifelse(is.na(qc_anno) & in_tissue, "None", qc_anno)),
           scran_qc_anno = factor(ifelse(as.logical(scran_discard), paste0("scran-", qc_anno), as.character(qc_anno)))
    )

pd|> count(in_tissue, qc_anno, scran_qc_anno)

spe$qc_anno <- pd$qc_anno
spe$scran_qc_anno <- pd$scran_qc_anno

levels(pd$scran_qc_anno)

## add to spe
table(spe$qc_anno)
table(spe$scran_anno)

qc_colors <- RColorBrewer::brewer.pal(name = "Set1", n=6)

qc_colors <- c(out_edge = "#E41A1C",
  fold = "#377EB8",
  out_tissue ="#4DAF4A",
  tail = "#984EA3",
  torn_edge ="#FF7F00",
  circle = "#FFFF33",
  None = "grey50")

qc_colors_dark <- c(out_edge = "#A51215",
               fold = "#26567E",
               out_tissue ="#327330",
               tail = "#67356E",
               torn_edge ="#B85C00",
               circle = "#CCCC00",
               None = "grey25")

names(qc_colors_dark) <- paste0("scran-", names(qc_colors_dark))


pdf(here(plot_dir, "spe_erc_QC_annotations.pdf"))
map(sample_order, 
    ~vis_clus(
        spe = spe,
        sampleid = .x,
        point_size = 1,
        clustervar = "scran_qc_anno",
        colors = c(qc_colors, qc_colors_dark)
    ))
dev.off()

