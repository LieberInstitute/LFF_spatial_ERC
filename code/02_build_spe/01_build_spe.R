
# cd /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/
suppressPackageStartupMessages(library("here"))
# remotes::install_github("drighelli/SpatialExperiment")
# remotes::install_github("LieberInstitute/spatialLIBD")
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Define some info for the samples
sample_info <- data.frame(
  sample_id = c(
    "V13Y24-343_A1",
    "V13Y24-343_B1",
    "V13Y24-343_C1",
    "V13Y24-343_D1",
    "V13Y24-344_A1",
    "V13Y24-344_B1",
    "V13Y24-344_C1",
    "V13Y24-344_D1",
    "V13Y24-342_A1",
    "V13Y24-342_B1",
    "V13Y24-342_C1",
    "V13Y24-342_D1",
    "V13Y24-340_A1",
    "V13Y24-340_B1",
    "V13Y24-340_C1",
    "V13Y24-340_D1",
    "V13B23-363_A1",
    "V13B23-363_B1",
    "V13B23-363_C1",
    "V13B23-363_D1",
    "V13B23-364_A1",
    "V13B23-364_B1",
    "V13B23-364_C1",
    "V13B23-364_D1",
    "V13B23-365_A1",
    "V13B23-365_B1",
    "V13B23-365_C1",
    "V13B23-365_D1",
    "V13B23-366_B1",
    "V13B23-366_C1",
    "V13B23-366_D1"
  )
)
sample_info$subject <- sample_info$sample_id
sample_info$sample_path <-
  file.path(
    here::here("processed-data", "01_spaceranger"),
    sample_info$sample_id,
    "outs"
  )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/raw-data/sample_info/Visium_HPC_Round1_20220113_Master_ADR.xlsx
## https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/raw-data/sample_info/Visium_HPC_Round2_20220223_Master_ADR.xlsx
donor_info <- data.frame(
  subject = c("V13Y24-343-A1","V13Y24-343-B1","V13Y24-343-C1","V13Y24-343-D1", "V13Y24-344-A1","V13Y24-344-B1","V13Y24-344-C1","V13Y24-344-D1","V13Y24-342-A1", "V13Y24-342-B1","V13Y24-342-C1","V13Y24-342-D1","V13Y24-340-A1","V13Y24-340-B1", "V13Y24-340-C1","V13Y24-340-D1", "V13B23-363_A1", "V13B23-363_B1", "V13B23-363_C1", "V13B23-363_D1", "V13B23-364_A1", "V13B23-364_B1", "V13B23-364_C1", "V13B23-364_D1", "V13B23-365_A1", "V13B23-365_B1", "V13B23-365_C1", "V13B23-365_D1", "V13B23-366_B1", "V13B23-366_C1", "V13B23-366_D1"),
  age = c(51.63, 51.45, 29.95, 59.86, 41.44, 60.83, 62.70, 46.53, 42.39, 48.75, 50.08, 61.34, 63.98, 54.88, 67.75, 31.31, 42.19, 55.88, 45.3, 51.11, 59.98, 43.67, 68.38, 58.19, 51.73, 50.2, 54.43, 48.69, 57.1, 60.84, 50.73),
  sex = c("M","F","F","M", "M", "M", "M", "F", "F", "F", "M", "M", "M", "F", "M", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F","M", "M", "F",),
  race = c("EA/CAUC", "AA", "AA", "EA/CAUC", "EA/CAUC", "AA", "AA", "AA", "AA", "AA", "EA/CAU", "AA", "EA/CAUC", "AA", "EA/CAUC", "EA/CAUC", "AA", "EA/CAUC", "AA", "EA/CAUC", "AA", "EA/CAUC", "AA", "AA", "EA/CAU", "AA", "EA/CAU", "AA", "EA/CAU", "EA/CAU", "AA"),
  diagnosis = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control"),
  rin = c(8.5, 8.6, 9.3, 7.4, 7.3, 8.7, 7.2, 8.5, 7.1, 9.2, 8.2, 8.5, 9.4, 8, 6.8, 6.6, 9, 8.5, 1, 8.6, 9.3, 5.2, 9, 8.7, 7.1, 8.3, 8.5, 5.2, 7.7, 7.4, 7), # fix rin for the donor with placeholder 1
  apoe = c("E2/E3", "E4/E4", "E2/E2", "E3/E4", "E2/E2", "E2/E3","E3/E4","E4/E4", "E3/E4", "E2/E3", "E4/E4", "E3/E4", "E3/ E4", "E4/E4", "E2/E2", "E2/E3", "E4/E4", "E3/E4", "E2/E2", "E2/E3", "E3/E4", "E2/E2", "E2/E3", "E4/E4", "E2/E3", "E3/E4", "E4/E4", "E2/E2", "E3/E4", "E2/E3", "E3/E4")
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE
)
Sys.time()

# [1] "2023-10-12 11:45:00 EDT"
# 2023-10-12 11:45:01.485023 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2023-10-12 11:52:25.451741 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2023-10-12 11:53:06.048223 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2023-10-12 11:53:14.299531 rtracklayer::import: reading the reference GTF file
# 2023-10-12 11:54:37.587021 adding gene information to the SPE object
# 2023-10-12 11:54:37.828386 adding information used by spatialLIBD
# [1] "2023-10-12 11:54:46 EDT"

## Add the study design info
add_design <- function(spe) {
  new_col <- merge(colData(spe), sample_info)
  ## Fix order
  new_col <- new_col[match(spe$key, new_col$key), ]
  stopifnot(identical(new_col$key, spe$key))
  rownames(new_col) <- rownames(colData(spe))
  colData(spe) <-
    new_col[, -which(colnames(new_col) == "sample_path")]
  return(spe)
}
spe <- add_design(spe)

# ## Read in cell counts and segmentation results
# segmentations_list <-
#     lapply(sample_info$sample_id, function(sampleid) {
#         file <-
#             here(
#                 "processed-data",
#                 "01_spaceranger_re-run",
#                 sampleid,
#                 "outs",
#                 "spatial",
#                 "tissue_spot_counts.csv"
#             )
#         if (!file.exists(file)) {
#             return(NULL)
#         }
#         x <- read.csv(file)
#         x$key <- paste0(x$barcode, "_", sampleid)
#         return(x)
#     })
#
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <-
#     Reduce(function(...) {
#         merge(..., all = TRUE)
#     }, segmentations_list[lengths(segmentations_list) > 0])

# ## Add the information
# segmentation_match <- match(spe$key, segmentations$key)
# segmentation_info <-
#     segmentations[segmentation_match, -which(
#         colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
#     )]
# colData(spe) <- cbind(colData(spe), segmentation_info)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 6675
length(no_expr) / nrow(spe) * 100
# [1] 18.23721
spe <- spe[-no_expr, ]


## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
  factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw <- spe

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

## Size in Gb
lobstr::obj_size(spe_raw)
# 4.86 GB

## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
# [1] 29926 122109
## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
  message("removing spots without counts for spe")
  spe <- spe[, -which(colSums(counts(spe)) == 0)]
  dim(spe)
}

# removing spots without counts for spe
# [1] 29926 122086

lobstr::obj_size(spe)
# 4.66 GB

saveRDS(spe, file.path(dir_rdata, "spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# > ## Reproducibility information
# print("Reproducibility information:")
# Sys.time()
# proc.time()
# options(width = 120)
# session_info()
# [1] "Reproducibility information:"
# [1] "2023-10-12 12:49:44 EDT"
#     user   system  elapsed
# 2284.341   30.228 4009.315
# 23-04-25 [2] Bioconductor
#  BiocFileCache            2.8.0     2023-04-25 [2] Bioconductor
#  BiocGenerics           * 0.46.0    2023-04-25 [2] Bioconductor
#  BiocIO                   1.10.0    2023-04-25 [2] Bioconductor
#  BiocManager              1.30.22   2023-08-08 [2] CRAN (R 4.3.1)
#  BiocNeighbors            1.18.0    2023-04-25 [2] Bioconductor
#  BiocParallel             1.34.2    2023-05-22 [2] Bioconductor
#  BiocSingular             1.16.0    2023-04-25 [2] Bioconductor
#  BiocVersion              3.17.1    2022-11-04 [2] Bioconductor
#  Biostrings               2.68.1    2023-05-16 [2] Bioconductor
#  bit                      4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
#  bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
#  bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
#  blob                     1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
#  bslib                    0.5.1     2023-08-11 [2] CRAN (R 4.3.1)
#  cachem                   1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
#  cli                      3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
#  codetools                0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
#  colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
#  config                   0.3.2     2023-08-30 [2] CRAN (R 4.3.1)
#  cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
#  crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
#  curl                     5.0.2     2023-08-14 [2] CRAN (R 4.3.1)
#  data.table               1.14.8    2023-02-17 [2] CRAN (R 4.3.1)
#  DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
#  dbplyr                   2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
#  DelayedArray             0.26.7    2023-07-28 [2] Bioconductor
#  DelayedMatrixStats       1.22.6    2023-08-28 [2] Bioconductor
#  digest                   0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
#  doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.3.1)
#  dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.3.1)
#  dplyr                    1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
#  dqrng                    0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
#  DropletUtils             1.20.0    2023-04-25 [2] Bioconductor
#  DT                       0.29      2023-08-29 [2] CRAN (R 4.3.1)
#  edgeR                    3.42.4    2023-05-31 [2] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.3.1)
#  ExperimentHub            2.8.1     2023-07-12 [2] Bioconductor
#  fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
#  fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
#  fields                   15.2      2023-08-17 [2] CRAN (R 4.3.1)
#  filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
#  foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.3.1)
#  generics                 0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
#  GenomeInfoDb           * 1.36.3    2023-09-07 [2] Bioconductor
#  GenomeInfoDbData         1.2.10    2023-07-20 [2] Bioconductor
#  GenomicAlignments        1.36.0    2023-04-25 [2] Bioconductor
#  GenomicRanges          * 1.52.0    2023-04-25 [2] Bioconductor
#  ggbeeswarm               0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
#  ggplot2                  3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
#  ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
#  glue                     1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
#  golem                    0.4.1     2023-06-05 [2] CRAN (R 4.3.1)
#  gridExtra                2.3       2017-09-09 [2] CRAN (R 4.3.1)
#  gtable                   0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
#  HDF5Array                1.28.1    2023-05-01 [2] Bioconductor
#  here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
#  htmltools                0.5.6     2023-08-10 [2] CRAN (R 4.3.1)
#  htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.3.1)
#  httpuv                   1.6.11    2023-05-11 [2] CRAN (R 4.3.1)
#  httr                     1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
#  interactiveDisplayBase   1.38.0    2023-04-25 [2] Bioconductor
#  IRanges                * 2.34.1    2023-06-22 [2] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
#  iterators                1.0.14    2022-02-05 [2] CRAN (R 4.3.1)
#  jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.3.1)
#  jsonlite                 1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
#  KEGGREST                 1.40.0    2023-04-25 [2] Bioconductor
#  later                    1.3.1     2023-05-02 [2] CRAN (R 4.3.1)
#  lattice                  0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.3.1)
#  lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
#  limma                    3.56.2    2023-06-04 [2] Bioconductor
#  lobstr                 * 1.1.2     2022-06-22 [2] CRAN (R 4.3.1)
#  locfit                   1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
#  magick                   2.7.5     2023-08-07 [2] CRAN (R 4.3.1)
#  magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
#  maps                     3.4.1     2022-10-30 [2] CRAN (R 4.3.1)
#  Matrix                   1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
#  MatrixGenerics         * 1.12.3    2023-07-30 [2] Bioconductor
#  matrixStats            * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
#  memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
#  mime                     0.12      2021-09-28 [2] CRAN (R 4.3.1)
#  munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
#  paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.3.1)
#  pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
#  pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
#  plotly                   4.10.2    2023-06-03 [2] CRAN (R 4.3.1)
#  png                      0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
#  prettyunits              1.1.1     2020-01-24 [2] CRAN (R 4.3.1)
#  promises                 1.2.1     2023-08-10 [2] CRAN (R 4.3.1)
#  purrr                    1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
#  R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
#  R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
#  R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
#  R6                       2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
#  rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
#  RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
#  Rcpp                     1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
#  RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
#  rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.3.1)
#  restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
#  rhdf5                    2.44.0    2023-04-25 [2] Bioconductor
#  rhdf5filters             1.12.1    2023-04-30 [2] Bioconductor
#  Rhdf5lib                 1.22.1    2023-09-10 [2] Bioconductor
#  rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
#  rlang                    1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
#  rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
#  Rsamtools                2.16.0    2023-04-25 [2] Bioconductor
#  RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
#  rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
#  rtracklayer            * 1.60.1    2023-08-15 [2] Bioconductor
#  S4Arrays                 1.0.6     2023-08-30 [2] Bioconductor
#  S4Vectors              * 0.38.1    2023-05-02 [2] Bioconductor
#  sass                     0.4.7     2023-07-15 [2] CRAN (R 4.3.1)
#  ScaledMatrix             1.8.1     2023-05-03 [2] Bioconductor
#  scales                   1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
#  scater                   1.28.0    2023-04-25 [2] Bioconductor
#  scuttle                  1.10.2    2023-08-03 [2] Bioconductor
#  sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
#  shiny                    1.7.5     2023-08-12 [2] CRAN (R 4.3.1)
#  shinyWidgets             0.8.0     2023-08-30 [2] CRAN (R 4.3.1)
#  SingleCellExperiment   * 1.22.0    2023-04-25 [2] Bioconductor
#  spam                     2.9-1     2022-08-07 [2] CRAN (R 4.3.1)
#  sparseMatrixStats        1.12.2    2023-07-02 [2] Bioconductor
#  SpatialExperiment      * 1.10.0    2023-04-25 [2] Bioconductor
#  spatialLIBD            * 1.12.0    2023-04-27 [2] Bioconductor
#  statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
#  SummarizedExperiment   * 1.30.2    2023-06-06 [2] Bioconductor
#  tibble                   3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
#  tidyr                    1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
#  tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
#  utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
#  vctrs                    0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
#  vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
#  viridis                  0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
#  viridisLite              0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
#  XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
#  xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.3.1)
#  XVector                  0.40.0    2023-04-25 [2] Bioconductor
#  yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
#  zlibbioc                 1.46.0    2023-04-25 [2] Bioconductor
#
#  [1] /users/hdivecha/R/4.3
#  [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
#  [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
#
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
