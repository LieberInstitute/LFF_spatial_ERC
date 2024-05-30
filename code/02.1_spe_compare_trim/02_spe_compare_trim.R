
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("rtracklayer")
library("lobstr")
library("tidyverse")
library("sessioninfo")

## Create output directories
plot_dir <- here::here("plots", "02.1_spe_compare_trim", "02_spe_comapre_trim")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
spe_trim <- readRDS(here("processed-data", "02.1_spe_compare_trim", "spe_trim.rds"))

pd <- as.data.frame(colData(spe_trim))

table(pd$sample, pd$trimmed)

pd |> filter(!trimmed) |> head()
pd |> filter(!trimmed) |> head()



trim_qc_wide <- pd |>
    select(sample, in_tissue, key, sum_umi, sum_gene, expr_chrM, expr_chrM_ratio, trimmed) |>
    mutate(trimmed = ifelse(trimmed, "trimmed", "untrimmed"),
           in_tissue = ifelse(in_tissue, "in_tissue", "out_tissue"),
           key = gsub("_untrimmed", "", key)) |>
    pivot_longer(!c(sample, in_tissue, key, trimmed), names_to = "metric", values_to = "value") |>
    pivot_wider(names_from = "trimmed", values_from = "value") |>
    mutate(diff = untrimmed - trimmed) |>
    separate(sample, into = c("slide", "section"), sep = "_", remove = FALSE)


trim_qc_wide |>
    group_by(metric) |>
    # summarise(as_tibble_row(quantile(diff)), median = median(diff))
    summarise(min_diff = min(diff),
              q25 = quantile(diff, .25),
              median_diff = median(diff),
              mean_diff = mean(diff),
              q25 = quantile(diff, .75),
              max_diff = max(diff))

# metric          min_diff       q25 median_diff mean_diff  max_diff
# 1 expr_chrM          -2     29         15        19.7       308     
# 2 expr_chrM_ratio    -0.05   0.00172    0.000632  0.000888    0.0455
# 3 sum_gene           -6     31.2       18        21.3       227     
# 4 sum_umi            -9    109         54        78.8      1718   

trim_qc_wide |>
    group_by(in_tissue, metric) |>
    # summarise(as_tibble_row(quantile(diff)), median = median(diff))
    summarise(min_diff = min(diff),
              q25 = quantile(diff, .25),
              median_diff = median(diff),
              mean_diff = mean(diff),
              q75 = quantile(diff, .75),
              max_diff = max(diff),
              sd_diff = sd(diff))

#   in_tissue  metric          min_diff       q25 median_diff mean_diff       q75  max_diff  sd_diff
# 1 in_tissue  expr_chrM        -1       8          19        22.9       32        308      19.9    
# 2 in_tissue  expr_chrM_ratio  -0.05   -0.000167    0.000619  0.000789   0.00156    0.0455  0.00196
# 3 in_tissue  sum_gene         -6      11          22        24.5       34        227      16.9    
# 4 in_tissue  sum_umi          -6      32          67        91.8      122       1718      91.1    
# 5 out_tissue expr_chrM        -2       0           1         1.54       2         50       2.30   
# 6 out_tissue expr_chrM_ratio  -0.0364 -0.00135     0.000874  0.00144    0.00375    0.0340  0.00440
# 7 out_tissue sum_gene         -5       1           3         3.49       5         52       3.92   
# 8 out_tissue sum_umi          -9       2           4         5.55       7        151       8.75      

trim_qc_wide |>
    group_by(in_tissue, metric) |>
    summarise(
              median_diff = median(diff),
              median_untr = median(untrimmed),
              median_precent_diff = 100*median_diff/median_untr,
              sd_diff = sd(diff),
              )

trim_qc_wide |>
    group_by(metric) |>
    count(diff < 0) 

## basic plots ##

qc_metrics_scatter <- trim_qc_wide |>
    ggplot(aes(trimmed, untrimmed)) +
    geom_point(aes(color = slide, shape = section), size = 1, alpha = 0.5) +
    geom_abline() +
    # facet_wrap(~metric, scales = "free")
    facet_wrap(in_tissue~metric, scales = "free", nrow = 2)

ggsave(qc_metrics_scatter, filename = here(plot_dir, "trim_qc_metrics_scatter.png"), width = 10)

qc_metrics_scatter_diff <- trim_qc_wide |>
    ggplot(aes(untrimmed, diff)) +
    geom_point(aes(color = slide, shape = section), size = 1, alpha = 0.5) +
    # facet_wrap(~metric, scales = "free")  +
    facet_wrap(in_tissue~metric, scales = "free", nrow = 2)+
    labs(y = "trimmed - untrimmed")

ggsave(qc_metrics_scatter_diff, filename = here(plot_dir, "trim_qc_metrics_scatter_diff.png"), width = 10)

qc_diff_boxplot <- trim_qc_wide |>
    ggplot(aes(sample, diff, color = slide)) +
    geom_boxplot() +
    # facet_wrap(~metric, scales = "free") +
    facet_wrap(in_tissue~metric, scales = "free", nrow = 2)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "trimmed - untrimmed")

ggsave(qc_diff_boxplot, filename = here(plot_dir, "trim_qc_diff_boxplot.png"), height = 8, width = 12)


qc_diff_boxplot_in_tissue <- trim_qc_wide |>
    ggplot(aes(in_tissue, diff)) +
    geom_boxplot() +
    facet_wrap(~metric, scales = "free") +
    # facet_wrap(in_tissue~metric, scales = "free", nrow = 2)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "trimmed - untrimmed")

ggsave(qc_diff_boxplot_in_tissue, filename = here(plot_dir, "trim_qc_diff_boxplot_in_tissue.png"), height = 8, width = 12)

# slurmjobs::job_single('02_spe_compare_trim', create_shell = TRUE, memory = '10G', command = "Rscript 02_spe_compare_trim.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()