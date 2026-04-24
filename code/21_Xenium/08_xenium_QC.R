## Louise Huuki-Myers, April 2026
## Build SummarizedExperiment Objetct for Xenium data
## Adapted from https://github.com/LieberInstitute/xenium_NAC/blob/54a6bd78ef1127d0e41f037bd6a080579dad9020/code/03_QC/01_qc.R

#### Set up ####

library("here")
library("SpatialExperiment")
library("qs2")

data_dir <- here("processed-data", "21_Xenium", "08_xenium_QC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), "- Loda data")
spe <- qs_read(here("processed-data", "21_Xenium", "07_xenium_build_spe","spe_xenium.qs2"))

spe
# dim: 541 473611

# Remove empty cells
message("Remove empty cells:", sum(colSums(counts(spe)) == 0)) #659
spe <- spe[, colSums(counts(spe)) > 0]           

#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_anyneg <- is_neg | is_neg2 | is_unassigned
is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel


