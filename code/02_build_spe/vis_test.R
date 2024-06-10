
library("spatialLIBD")
library("here")

plot_dir <- here("plots", "02_build_spe", "vis_test")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

spe <- readRDS(here("processed-data", "02_build_spe","spe_raw.rds"))

table(spe$sample_id, spe$in_tissue)
#        FALSE TRUE
# Br1039  1073 3919
# Br1289  1329 3663
# Br1556  1458 3534
# Br1691   806 4186
# Br1706  1238 3754
# Br2305  1910 3082
# Br2582  2596 2396
# Br3974  2264 2728
# Br5161   292 4700
# Br5212   315 4677
# ...


colData(spe)

vis_grid_clus(
    spe = spe,
    clustervar = "in_tissue",
    pdf = here::here(plot_dir, "spe_erc_grid-in_tissue.pdf"),
    sort_clust = FALSE,
    point_size = 1
    # colors = c("in" = "grey90", "out" = "orange")
)
# With sort_clust = TRUE
# Error in `[[<-`(`*tmp*`, clustervar, value = c(`FALSE` = 2L, `TRUE` = 1L,  : 
#                                                    122746 elements in value to replace 154752 elements

pdf(here::here(plot_dir, "spe_erc_Br6085-in_tissue.pdf"))
vis_clus(
    spe = spe,
    sample_id = "Br6085",
    clustervar = "in_tissue",
    # clustervar = "ManualAnnotation",
    # sort_clust = FALSE,
    point_size = 1
    # colors = c("in" = "grey90", "out" = "orange")
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.3 (2024-02-29)
# os       macOS Sonoma 14.5
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-06-10
# rstudio  2024.04.1+748 Chocolate Cosmos (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version     date (UTC) lib source
# abind                    1.4-5       2016-07-21 [1] CRAN (R 4.3.0)
# AnnotationDbi            1.64.1      2023-11-03 [1] Bioconductor
# AnnotationHub            3.10.1      2024-04-05 [1] Bioconductor 3.18 (R 4.3.3)
# attempt                  0.3.1       2020-05-03 [1] CRAN (R 4.3.0)
# beachmat                 2.18.1      2024-02-14 [1] Bioconductor 3.18 (R 4.3.2)
# beeswarm                 0.4.0       2021-06-01 [1] CRAN (R 4.3.0)
# benchmarkme              1.0.8       2022-06-12 [1] CRAN (R 4.3.0)
# benchmarkmeData          1.0.4       2020-04-23 [1] CRAN (R 4.3.0)
# Biobase                * 2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache            2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.3)
# BiocGenerics           * 0.48.1      2023-11-01 [1] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [1] Bioconductor
# BiocManager              1.30.23     2024-05-04 [1] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0      2023-10-24 [1] Bioconductor
# BiocSingular             1.18.0      2023-10-24 [1] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [1] Bioconductor
# Biostrings               2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.3)
# bit                      4.0.5       2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5       2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7       2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4       2023-03-17 [1] CRAN (R 4.3.0)
# bslib                    0.7.0       2024-03-29 [1] CRAN (R 4.3.2)
# cachem                   1.1.0       2024-05-16 [1] CRAN (R 4.3.3)
# cli                      3.6.2       2023-12-11 [1] CRAN (R 4.3.0)
# codetools                0.2-20      2024-03-31 [1] CRAN (R 4.3.2)
# colorout               * 1.3-0.2     2024-03-21 [1] Github (jalvesaq/colorout@c6113a2)
# colorspace               2.1-0       2023-01-23 [1] CRAN (R 4.3.0)
# config                   0.3.2       2023-08-30 [1] CRAN (R 4.3.0)
# cowplot                  1.1.3       2024-01-22 [1] CRAN (R 4.3.2)
# crayon                   1.5.2       2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.1       2024-03-01 [1] CRAN (R 4.3.2)
# data.table               1.15.4      2024-03-30 [1] CRAN (R 4.3.2)
# DBI                      1.2.3       2024-06-02 [1] CRAN (R 4.3.3)
# dbplyr                   2.5.0       2024-03-19 [1] CRAN (R 4.3.2)
# DelayedArray             0.28.0      2023-10-24 [1] Bioconductor
# DelayedMatrixStats       1.24.0      2023-10-24 [1] Bioconductor
# devtools               * 2.4.5       2022-10-11 [1] CRAN (R 4.3.0)
# digest                   0.6.35      2024-03-11 [1] CRAN (R 4.3.2)
# doParallel               1.0.17      2022-02-07 [1] CRAN (R 4.3.0)
# dotCall64                1.1-1       2023-11-28 [1] CRAN (R 4.3.0)
# dplyr                    1.1.4       2023-11-17 [1] CRAN (R 4.3.0)
# DT                       0.33        2024-04-04 [1] CRAN (R 4.3.2)
# edgeR                    4.0.16      2024-02-18 [1] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [1] CRAN (R 4.3.0)
# ExperimentHub            2.10.0      2023-10-24 [1] Bioconductor
# fansi                    1.0.6       2023-12-08 [1] CRAN (R 4.3.0)
# farver                   2.1.2       2024-05-13 [1] CRAN (R 4.3.3)
# fastmap                  1.2.0       2024-05-15 [1] CRAN (R 4.3.3)
# fields                   15.2        2023-08-17 [1] CRAN (R 4.3.0)
# filelock                 1.0.3       2023-12-11 [1] CRAN (R 4.3.0)
# foreach                  1.5.2       2022-02-02 [1] CRAN (R 4.3.0)
# fs                       1.6.4       2024-04-25 [1] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.3)
# GenomeInfoDbData         1.2.11      2024-03-21 [1] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [1] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [1] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [1] CRAN (R 4.3.0)
# ggplot2                  3.5.1       2024-04-23 [1] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [1] CRAN (R 4.3.0)
# glue                     1.7.0       2024-01-09 [1] CRAN (R 4.3.0)
# golem                    0.4.1       2023-06-05 [1] CRAN (R 4.3.0)
# gridExtra                2.3         2017-09-09 [1] CRAN (R 4.3.0)
# gtable                   0.3.5       2024-04-22 [1] CRAN (R 4.3.2)
# here                   * 1.0.1       2020-12-13 [1] CRAN (R 4.3.0)
# htmltools                0.5.8.1     2024-04-04 [1] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [1] CRAN (R 4.3.0)
# httpuv                   1.6.15      2024-03-26 [1] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [1] CRAN (R 4.3.0)
# interactiveDisplayBase   1.40.0      2023-10-24 [1] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [1] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [1] CRAN (R 4.3.0)
# iterators                1.0.14      2022-02-05 [1] CRAN (R 4.3.0)
# jquerylib                0.1.4       2021-04-26 [1] CRAN (R 4.3.0)
# jsonlite                 1.8.8       2023-12-04 [1] CRAN (R 4.3.0)
# KEGGREST                 1.42.0      2023-10-24 [1] Bioconductor
# labeling                 0.4.3       2023-08-29 [1] CRAN (R 4.3.0)
# later                    1.3.2       2023-12-06 [1] CRAN (R 4.3.0)
# lattice                  0.22-6      2024-03-20 [1] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.4       2023-11-07 [1] CRAN (R 4.3.0)
# limma                    3.58.1      2023-10-31 [1] Bioconductor
# locfit                   1.5-9.9     2024-03-01 [1] CRAN (R 4.3.2)
# magick                   2.8.3       2024-02-18 [1] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [1] CRAN (R 4.3.0)
# maps                     3.4.2       2023-12-15 [1] CRAN (R 4.3.0)
# Matrix                   1.6-5       2024-01-11 [1] CRAN (R 4.3.3)
# MatrixGenerics         * 1.14.0      2023-10-24 [1] Bioconductor
# matrixStats            * 1.3.0       2024-04-11 [1] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12        2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                   0.1.1.1     2018-05-18 [1] CRAN (R 4.3.0)
# munsell                  0.5.1       2024-04-01 [1] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [1] CRAN (R 4.3.0)
# pillar                   1.9.0       2023-03-22 [1] CRAN (R 4.3.0)
# pkgbuild                 1.4.4       2024-03-17 [1] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [1] CRAN (R 4.3.0)
# pkgload                  1.3.4       2024-01-16 [1] CRAN (R 4.3.0)
# plotly                   4.10.4      2024-01-13 [1] CRAN (R 4.3.0)
# png                      0.1-8       2022-11-29 [1] CRAN (R 4.3.0)
# profvis                  0.3.8       2023-05-02 [1] CRAN (R 4.3.0)
# promises                 1.3.0       2024-04-05 [1] CRAN (R 4.3.2)
# purrr                    1.0.2       2023-08-10 [1] CRAN (R 4.3.0)
# R6                       2.5.1       2021-08-19 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3       2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3       2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.12      2024-01-09 [1] CRAN (R 4.3.0)
# RCurl                    1.98-1.14   2024-01-09 [1] CRAN (R 4.3.0)
# rematch2                 2.1.2       2020-05-01 [1] CRAN (R 4.3.0)
# remotes                  2.5.0       2024-03-17 [1] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [1] CRAN (R 4.3.0)
# rjson                    0.2.21      2022-01-09 [1] CRAN (R 4.3.0)
# rlang                    1.1.4       2024-06-04 [1] CRAN (R 4.3.3)
# rprojroot                2.0.4       2023-11-05 [1] CRAN (R 4.3.0)
# Rsamtools                2.18.0      2023-10-24 [1] Bioconductor
# RSQLite                  2.3.7       2024-05-27 [1] CRAN (R 4.3.3)
# rstudioapi               0.16.0      2024-03-24 [1] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [1] CRAN (R 4.3.0)
# rtracklayer              1.62.0      2023-10-24 [1] Bioconductor
# S4Arrays                 1.2.1       2024-03-06 [1] Bioconductor 3.18 (R 4.3.3)
# S4Vectors              * 0.40.2      2023-11-23 [1] Bioconductor
# sass                     0.4.9       2024-03-15 [1] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [1] Bioconductor
# scales                   1.3.0       2023-11-28 [1] CRAN (R 4.3.0)
# scater                   1.30.1      2023-12-06 [1] Bioconductor
# scuttle                  1.12.0      2023-10-24 [1] Bioconductor
# sessioninfo              1.2.2       2021-12-06 [1] CRAN (R 4.3.0)
# shiny                    1.8.1.1     2024-04-02 [1] CRAN (R 4.3.2)
# shinyWidgets             0.8.6       2024-04-24 [1] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [1] Bioconductor
# spam                     2.10-0      2023-10-23 [1] CRAN (R 4.3.0)
# SparseArray              1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [1] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [1] Bioconductor
# spatialLIBD            * 1.17.2      2024-06-06 [1] Github (LieberInstitute/spatialLIBD@e0947a2)
# statmod                  1.5.0       2023-01-06 [1] CRAN (R 4.3.0)
# stringi                  1.8.4       2024-05-06 [1] CRAN (R 4.3.2)
# stringr                  1.5.1       2023-11-14 [1] CRAN (R 4.3.0)
# SummarizedExperiment   * 1.32.0      2023-10-24 [1] Bioconductor
# tibble                   3.2.1       2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                    1.3.1       2024-01-24 [1] CRAN (R 4.3.2)
# tidyselect               1.2.1       2024-03-11 [1] CRAN (R 4.3.2)
# urlchecker               1.0.1       2021-11-30 [1] CRAN (R 4.3.0)
# usethis                * 2.2.3       2024-02-19 [1] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [1] CRAN (R 4.3.0)
# vctrs                    0.6.5       2023-12-01 [1] CRAN (R 4.3.0)
# vipor                    0.4.7       2023-12-18 [1] CRAN (R 4.3.0)
# viridis                  0.6.5       2024-01-29 [1] CRAN (R 4.3.2)
# viridisLite              0.4.2       2023-05-02 [1] CRAN (R 4.3.0)
# withr                    3.0.0       2024-01-16 [1] CRAN (R 4.3.0)
# XML                      3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.2)
# xtable                   1.8-4       2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0      2023-10-24 [1] Bioconductor
# yaml                     2.3.8       2023-12-11 [1] CRAN (R 4.3.0)
# zlibbioc                 1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.3)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
