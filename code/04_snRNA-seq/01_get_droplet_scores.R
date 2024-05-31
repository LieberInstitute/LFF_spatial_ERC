# Collaborated with Louise Huuki 
# Based on
# https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R
# list.files("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger")
# qrsh -l mem_free=75G,h_vmem=75G

library("SingleCellExperiment")
library("DropletUtils")
# library("BiocParallel")
library("scuttle")
library("tidyverse")
library("here")
library("sessioninfo")

## get sample i
args <- commandArgs(trailingOnly = TRUE)
sample <- args[[1]]
sample_path <- here("processed-data", "07_cellranger", sample, "outs", "raw_feature_bc_matrix") 
stopifnot(file.exists(sample_path))

message(Sys.time(), " Reading data from ", sample_path)

#### Load & Subset raw data ####
sce <- read10xCounts(sample_path, col.names=TRUE) 

message("ncol:", ncol(sce))

#### Run barcodeRanks to find knee ####
message(Sys.time(), "Running barcode ranks.")

bcRanks <- barcodeRanks(sce, fit.bounds = c(10, 1e3))

knee_lower <- metadata(bcRanks)$knee + 100
message(
  "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee_lower =", knee_lower
)

#### Run emptyDrops w/ knee + 100 ####
set.seed(100)
message(Sys.time(), "Starting emptyDrops")
e.out <- DropletUtils::emptyDrops(
  sce,
  niters = 30000,
  lower = knee_lower
  # ,
  # BPPARAM = BiocParallel::MulticoreParam(4)
)
message(Sys.time(), "Done - saving data")

save(e.out, file = here("processed-data", "04_snRNA-seq", "01_get_droplet_scores", paste0("droplet_scores_", sample, ".Rdata")))

#### QC Plots ####
message("QC check")
FDR_cutoff <- 0.001
addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))

n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
message(n_cell_anno)

my_theme <- theme_bw() +
  theme(text = element_text(size = 15))

droplet_elbow_plot <- as.data.frame(bcRanks) %>%
  add_column(FDR = e.out$FDR) %>%
  ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
  annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Second Knee", vjust = -1, color = "gray") +
  geom_hline(yintercept = knee_lower, linetype = "dashed") +
  annotate("text", x = 10, y = knee_lower, label = "Knee est 'lower'", vjust = -0.5) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Barcode Rank",
    y = "Total UMIs",
    title = paste("Sample", sample),
    subtitle = n_cell_anno,
    color = paste("FDR <", FDR_cutoff)
  ) +
  my_theme +
  theme(legend.position = "bottom")

# droplet_scatter_plot <- as.data.frame(e) %>%
#   ggplot(aes(x = Total, y = -LogProb, color = FDR < FDR_cutoff)) +
#   geom_point(alpha = 0.5, size = 1) +
#   labs(x = "Total UMIs", y = "-Log Probability",
#        color = paste("FDR <", FDR_cutoff)) +
#   my_theme+
#   theme(legend.position = "bottom")
# # print(droplet_elbow_plot/droplet_scatter_plot)
# ggsave(droplet_elbow_plot/droplet_scatter_plot, filename = here("plots","03_build_sce", "droplet_qc_png",paste0("droplet_qc_",sample,".png")))

ggsave(droplet_elbow_plot, filename = here("plots", "04_snRNA-seq", "01_get_droplet_scores", paste0("droplet_qc_", sample, ".png")))


# sgejobs::job_single('01_get_droplet_scores', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 01_get_droplet_scores.R Br1092")
# sample_list <- list.files("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/07_cellranger")
# sgejobs::job_loop(loops = list(sample = sample_list), name = "01_get_droplet_scores_loop", 
#          create_shell = TRUE, queue = "bluejay", memory = "50G")



## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.3 Patched (2023-04-07 r84211)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-06-13
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# beachmat               2.14.2    2023-04-07 [2] Bioconductor
# Biobase              * 2.58.0    2022-11-01 [2] Bioconductor
# BiocGenerics         * 0.44.0    2022-11-01 [2] Bioconductor
# BiocParallel           1.32.6    2023-03-17 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils         * 1.18.1    2022-11-22 [2] Bioconductor
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# ggplot2              * 3.4.2     2023-04-03 [2] CRAN (R 4.2.3)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.2.3)
# HDF5Array              1.26.0    2022-11-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.2.3)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.2.3)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# purrr                * 1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.2.2)
# rhdf5                  2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters           1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib               1.20.0    2022-11-01 [2] Bioconductor
# rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scuttle              * 1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.2.2)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.2.2)
# tzdb                   0.3.0     2022-03-28 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.2     2023-04-19 [1] CRAN (R 4.2.3)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# XVector                0.38.0    2022-11-01 [2] Bioconductor
# zlibbioc               1.44.0    2022-11-01 [2] Bioconductor
# 
# [1] /users/bsimbiat/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
