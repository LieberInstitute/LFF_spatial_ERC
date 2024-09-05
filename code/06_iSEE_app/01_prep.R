library("SingleCellExperiment")
library("here")
library("lobstr")
library("sessioninfo")

## Load the sce data
load(here("processed-data", "sce_objects", "sce_Habenula_Pilot.Rdata"), verbose = TRUE)

## Drop data we don't need for iSEE
assayNames(sce)
# [1] "counts"                      "binomial_deviance_residuals" "logcounts"
assays(sce)$counts <- NULL
assays(sce)$binomial_deviance_residuals <- NULL

## Drop metadata we don't need
metadata(sce) <- list()
sce$path <- NULL
sce$total <- NULL

# sourcing official color palette
source(
    file = here("code", "99_paper_figs", "source_colors.R"),
    echo = TRUE
)

## Check final size
lobstr::obj_size(sce)
# 988.31 MB

saveRDS(sce, file = here("code", "16_iSEE", "02_snRNAseq", "sce.rds"))
saveRDS(sn_colors, file = here("code", "16_iSEE", "02_snRNAseq", "sn_colors.rds"))
saveRDS(bulk_colors, file = here("code", "16_iSEE", "02_snRNAseq", "bulk_colors.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 (2023-06-16)
#  os       macOS Sonoma 14.2.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2024-01-26
#  rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
#  pandoc   3.1.5 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
#  Biobase              * 2.60.0    2023-04-25 [1] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [1] Bioconductor
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  brio                   1.1.4     2023-12-10 [1] CRAN (R 4.3.1)
#  cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
#  colorout               1.3-0     2023-09-28 [1] Github (jalvesaq/colorout@8384882)
#  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  data.table             1.14.10   2023-12-08 [1] CRAN (R 4.3.1)
#  DelayedArray           0.26.7    2023-07-30 [1] Bioconductor
#  devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
#  digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
#  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.4    2023-10-02 [1] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-05-06 [1] Bioconductor
#  GenomicRanges        * 1.52.1    2023-10-08 [1] Bioconductor
#  glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
#  htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.1)
#  httpuv                 1.6.13    2023-12-06 [1] CRAN (R 4.3.1)
#  IRanges              * 2.34.1    2023-07-02 [1] Bioconductor
#  later                  1.3.2     2023-12-06 [1] CRAN (R 4.3.1)
#  lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
#  lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
#  lobstr                 1.1.2     2022-06-22 [1] CRAN (R 4.3.0)
#  lubridate              1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.3    2023-07-30 [1] Bioconductor
#  matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  pkgbuild               1.4.3     2023-12-10 [1] CRAN (R 4.3.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                1.3.4     2024-01-16 [1] CRAN (R 4.3.1)
#  prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
#  profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                 1.0.2     2023-08-31 [1] CRAN (R 4.3.0)
#  purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
#  RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
#  remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
#  rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
#  rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
#  rsthemes               0.4.0     2023-05-06 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
#  S4Arrays               1.0.6     2023-09-10 [1] Bioconductor
#  S4Vectors            * 0.38.2    2023-09-22 [1] Bioconductor
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  shiny                  1.8.0     2023-11-17 [1] CRAN (R 4.3.1)
#  SingleCellExperiment * 1.22.0    2023-04-25 [1] Bioconductor
#  stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
#  stringr                1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
#  SummarizedExperiment * 1.30.2    2023-06-11 [1] Bioconductor
#  suncalc                0.5.1     2022-09-29 [1] CRAN (R 4.3.0)
#  testthat             * 3.2.1     2023-12-02 [1] CRAN (R 4.3.1)
#  timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
#  urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
#  usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
#  vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [1] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
