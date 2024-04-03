library("jsonlite")
library("here")
library("stringr")
library("sessioninfo")

params_part1 <- read.csv(here("code", "01_spaceranger", "spaceranger_parameters_1v-16v.txt"), header = FALSE)
params_part2 <- read.csv(here("code", "01_spaceranger", "spaceranger_parameters_17v-31v.txt"), header = FALSE)
params <- rbind(params_part1, params_part2)
colnames(params) <- c("slide_capture_area", "slide", "capture_area", "input_image", "json", "fastq")

for(i in seq_len(nrow(params))) {
	message("We are editing json file ", params$json[i])
	json_raw <- readLines(params$json[i])
	print(stringr::str_sub(json_raw, start= -100))
	json_raw <- gsub('\\"serialNumber\\"\\:\\"\\",\\"area\\"\\:\\"\\"', paste0('\\"serialNumber\\"\\:\\"', params$slide[i], '\\",\\"area\\"\\:\\"', params$capture_area[i], '\\"'), json_raw)
	print(stringr::str_sub(json_raw, start= -110))
	writeLines(json_raw, params$json[i])
	
	json <- read_json(params$json[i])
	stopifnot(json$slideNumber == params$slide[i])
	stopifnot(json$area == params$capture_area[i])
}

# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-343-A1.json
# [1] "umber\":\"V13Y24-343\",\"area\":\"A1\",\"checksum\":\"acc29ab307ea86a3cc709ca9f2d55fe9\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-343\",\"area\":\"A1\",\"checksum\":\"acc29ab307ea86a3cc709ca9f2d55fe9\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-343-B1.json
# [1] "umber\":\"V13Y24-343\",\"area\":\"B1\",\"checksum\":\"86ef89e2065c4b1f11c04c8186bbfbed\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-343\",\"area\":\"B1\",\"checksum\":\"86ef89e2065c4b1f11c04c8186bbfbed\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-343-C1.json
# [1] "umber\":\"V13Y24-343\",\"area\":\"C1\",\"checksum\":\"d4d6cb88072b09e6ce7c64ab0459f98f\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-343\",\"area\":\"C1\",\"checksum\":\"d4d6cb88072b09e6ce7c64ab0459f98f\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-343-D1.json
# [1] "umber\":\"V13Y24-343\",\"area\":\"D1\",\"checksum\":\"8091eb85070d8c9a0a67fc30dd6c0890\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-343\",\"area\":\"D1\",\"checksum\":\"8091eb85070d8c9a0a67fc30dd6c0890\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-344-A1.json
# [1] "umber\":\"V13Y24-344\",\"area\":\"A1\",\"checksum\":\"d0fce7e4057108325207b6a0d6e6caeb\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-344\",\"area\":\"A1\",\"checksum\":\"d0fce7e4057108325207b6a0d6e6caeb\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-344-B1.json
# [1] "umber\":\"V13Y24-344\",\"area\":\"B1\",\"checksum\":\"fed58aa507a2f2afd0bd2e8a8c79506d\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-344\",\"area\":\"B1\",\"checksum\":\"fed58aa507a2f2afd0bd2e8a8c79506d\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-344-C1.json
# [1] "umber\":\"V13Y24-344\",\"area\":\"C1\",\"checksum\":\"d6eddbd307fd4b86715d8fb43d50b681\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-344\",\"area\":\"C1\",\"checksum\":\"d6eddbd307fd4b86715d8fb43d50b681\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-344-D1.json
# [1] "umber\":\"V13Y24-344\",\"area\":\"D1\",\"checksum\":\"c2f583d30749bdef0c1975bf2205963a\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-344\",\"area\":\"D1\",\"checksum\":\"c2f583d30749bdef0c1975bf2205963a\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-342-A1.json
# [1] "umber\":\"V13Y24-342\",\"area\":\"A1\",\"checksum\":\"f4a3dcea07e15d2b9e59320a78e1a096\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-342\",\"area\":\"A1\",\"checksum\":\"f4a3dcea07e15d2b9e59320a78e1a096\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-342-B1.json
# [1] "umber\":\"V13Y24-342\",\"area\":\"B1\",\"checksum\":\"e091a0e20726172264435c320047e901\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-342\",\"area\":\"B1\",\"checksum\":\"e091a0e20726172264435c320047e901\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-342-C1.json
# [1] "umber\":\"V13Y24-342\",\"area\":\"C1\",\"checksum\":\"96d2f8e907e8c2dadf26a354dd916588\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-342\",\"area\":\"C1\",\"checksum\":\"96d2f8e907e8c2dadf26a354dd916588\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-342-D1.json
# [1] "umber\":\"V13Y24-342\",\"area\":\"D1\",\"checksum\":\"caf6d6e107cbc4c934536fdc2bdc567d\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-342\",\"area\":\"D1\",\"checksum\":\"caf6d6e107cbc4c934536fdc2bdc567d\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-340-A1.json
# [1] "umber\":\"V13Y24-340\",\"area\":\"A1\",\"checksum\":\"c28e82810871f766bb383863d4cc813e\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-340\",\"area\":\"A1\",\"checksum\":\"c28e82810871f766bb383863d4cc813e\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-340-B1.json
# [1] "umber\":\"V13Y24-340\",\"area\":\"B1\",\"checksum\":\"94af8981088c4d82e97c0d481fe9db5b\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-340\",\"area\":\"B1\",\"checksum\":\"94af8981088c4d82e97c0d481fe9db5b\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-340-C1.json
# [1] "umber\":\"V13Y24-340\",\"area\":\"C1\",\"checksum\":\"c6a982f1075b5183543aa365fe86c856\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-340\",\"area\":\"C1\",\"checksum\":\"c6a982f1075b5183543aa365fe86c856\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13Y24-340-D1.json
# [1] "umber\":\"V13Y24-340\",\"area\":\"D1\",\"checksum\":\"e28ae6c807c3ec336058931292f615b8\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13Y24-340\",\"area\":\"D1\",\"checksum\":\"e28ae6c807c3ec336058931292f615b8\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-363-A1.json
# [1] "umber\":\"V13B23-363\",\"area\":\"A1\",\"checksum\":\"c48dc1e50708d1b78af9f9b85cacfdf5\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-363\",\"area\":\"A1\",\"checksum\":\"c48dc1e50708d1b78af9f9b85cacfdf5\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-363-B1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"f0a38ae607c16a7af817a1031274bf42\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-363\",\"area\":\"B1\",\"checksum\":\"f0a38ae607c16a7af817a1031274bf42\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-363-C1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"8aaea1e607b9799a8c068fb7135e656b\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-363\",\"area\":\"C1\",\"checksum\":\"8aaea1e607b9799a8c068fb7135e656b\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-363-D1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"e6d088e707a0e84da9f518cde33f0b9d\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-363\",\"area\":\"D1\",\"checksum\":\"e6d088e707a0e84da9f518cde33f0b9d\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-364-A1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"eecff488085fcf9b7b92cce77cac91af\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-364\",\"area\":\"A1\",\"checksum\":\"eecff488085fcf9b7b92cce77cac91af\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-364-B1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"caa0d084081121daf1218825a4f1f804\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-364\",\"area\":\"B1\",\"checksum\":\"caa0d084081121daf1218825a4f1f804\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-364-C1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"cec9b8870844783a21686d6da7fa102f\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-364\",\"area\":\"C1\",\"checksum\":\"cec9b8870844783a21686d6da7fa102f\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-364-D1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"caafef830882d8af600a3e1bffad8b8f\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-364\",\"area\":\"D1\",\"checksum\":\"caafef830882d8af600a3e1bffad8b8f\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-365-A1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"d6d391ba08f67f18e38c019ffe045c0e\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-365\",\"area\":\"A1\",\"checksum\":\"d6d391ba08f67f18e38c019ffe045c0e\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-365-B1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"fea6cfb70860b44a5a2e2d20370378c2\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-365\",\"area\":\"B1\",\"checksum\":\"fea6cfb70860b44a5a2e2d20370378c2\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-365-C1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"aefea0a908ed1d587e7d55225cd1d39e\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-365\",\"area\":\"C1\",\"checksum\":\"aefea0a908ed1d587e7d55225cd1d39e\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-365-D1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"c88389b408469cb684ed06f0a8298185\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-365\",\"area\":\"D1\",\"checksum\":\"c88389b408469cb684ed06f0a8298185\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-366-B1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"90b8cb900848234956c00f3cd61f2696\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-366\",\"area\":\"B1\",\"checksum\":\"90b8cb900848234956c00f3cd61f2696\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-366-C1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"94a495cc07cf61f155c15a025513c600\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-366\",\"area\":\"C1\",\"checksum\":\"94a495cc07cf61f155c15a025513c600\",\"removeImagePages\":[]}"
# We are editing json file /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/Images/loupe-alignment/V13B23-366-D1.json
# [1] "1]],\"serialNumber\":\"\",\"area\":\"\",\"checksum\":\"bab389d3053426ce30663ed4dae07ee8\",\"removeImagePages\":[]}"
# [1] "],\"serialNumber\":\"V13B23-366\",\"area\":\"D1\",\"checksum\":\"bab389d3053426ce30663ed4dae07ee8\",\"removeImagePages\":[]}"

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 Patched (2024-02-08 r85876)
#  os       Rocky Linux 9.2 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2024-04-03
#  pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────
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
#  stringi       1.8.2   2023-11-23 [1] CRAN (R 4.3.2)
#  stringr     * 1.5.1   2023-11-14 [2] CRAN (R 4.3.2)
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
# ─────────────────────────────────────────────────────────────────────────────────
