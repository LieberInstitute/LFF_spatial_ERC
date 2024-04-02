library("jsonlite")
library("here")
library("sessioninfo")

params <- read.csv(here("code", "01_spaceranger", "spaceranger_parameters.txt"), header = FALSE)
colnames(params) <- c("slide_capture_area", "slide", "capture_area", "input_image", "json", "fastq")

for(i in seq_len(nrow(params))) {
	message("We are editing json file ", params$json[i])
	json <- read_json(params$json[i])
	json$serialNumber <- params$slide[i]
	json$area <- params$capture_area[i]
	write_json(json, path = params$json[i])
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 Patched (2024-02-08 r85876)
#  os       Rocky Linux 9.2 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2024-04-02
#  pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date (UTC) lib source
#  cli           3.6.2   2023-12-11 [2] CRAN (R 4.3.2)
#  colorout    * 1.3-0.1 2023-11-29 [1] Github (jalvesaq/colorout@deda341)
#  colorspace    2.1-0   2023-01-23 [2] CRAN (R 4.3.2)
#  digest        0.6.34  2024-01-11 [2] CRAN (R 4.3.2)
#  dplyr         1.1.4   2023-11-17 [1] CRAN (R 4.3.2)
#  fansi         1.0.6   2023-12-08 [2] CRAN (R 4.3.2)
#  fastmap       1.1.1   2023-02-24 [2] CRAN (R 4.3.2)
#  generics      0.1.3   2022-07-05 [2] CRAN (R 4.3.2)
#  ggplot2       3.4.4   2023-10-12 [2] CRAN (R 4.3.2)
#  glue          1.7.0   2024-01-09 [2] CRAN (R 4.3.2)
#  gtable        0.3.4   2023-08-21 [2] CRAN (R 4.3.2)
#  here        * 1.0.1   2020-12-13 [2] CRAN (R 4.3.2)
#  htmltools     0.5.7   2023-11-03 [2] CRAN (R 4.3.2)
#  htmlwidgets   1.6.4   2023-12-06 [2] CRAN (R 4.3.2)
#  httpuv        1.6.14  2024-01-26 [2] CRAN (R 4.3.2)
#  jsonlite    * 1.8.8   2023-12-04 [2] CRAN (R 4.3.2)
#  later         1.3.2   2023-12-06 [2] CRAN (R 4.3.2)
#  lattice       0.22-5  2023-10-24 [3] CRAN (R 4.3.2)
#  lifecycle     1.0.4   2023-11-07 [2] CRAN (R 4.3.2)
#  magrittr      2.0.3   2022-03-30 [2] CRAN (R 4.3.2)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.3.2)
#  pillar        1.9.0   2023-03-22 [2] CRAN (R 4.3.2)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.3.2)
#  png           0.1-8   2022-11-29 [2] CRAN (R 4.3.2)
#  promises      1.2.1   2023-08-10 [2] CRAN (R 4.3.2)
#  R6            2.5.1   2021-08-19 [2] CRAN (R 4.3.2)
#  Rcpp          1.0.12  2024-01-09 [2] CRAN (R 4.3.2)
#  rlang         1.1.3   2024-01-10 [2] CRAN (R 4.3.2)
#  rmote         0.3.4   2023-11-29 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot     2.0.4   2023-11-05 [2] CRAN (R 4.3.2)
#  scales        1.3.0   2023-11-28 [2] CRAN (R 4.3.2)
#  servr         0.27    2023-05-02 [1] CRAN (R 4.3.2)
#  sessioninfo * 1.2.2   2021-12-06 [2] CRAN (R 4.3.2)
#  tibble        3.2.1   2023-03-20 [2] CRAN (R 4.3.2)
#  tidyselect    1.2.0   2022-10-10 [2] CRAN (R 4.3.2)
#  utf8          1.2.4   2023-10-22 [2] CRAN (R 4.3.2)
#  vctrs         0.6.5   2023-12-01 [2] CRAN (R 4.3.2)
#  xfun          0.42    2024-02-08 [2] CRAN (R 4.3.2)
# 
#  [1] /users/lcollado/R/4.3.x
#  [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
#  [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
