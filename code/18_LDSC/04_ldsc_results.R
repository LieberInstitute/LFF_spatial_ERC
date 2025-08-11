## Louise Huuki-Myers, Aug 2025
## Compile LDSC results
## based on https://github.com/LieberInstitute/spatialdACC/blob/main/code/17-2_LDSC_enrichment/spatial/ldsc_results2.R

#### Set Up ####

library("tidyverse")
library("getopt")
library("SingleCellExperiment")
library("here")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type",
      "mode", "m", "1", "character", "specificity or enrichment"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)