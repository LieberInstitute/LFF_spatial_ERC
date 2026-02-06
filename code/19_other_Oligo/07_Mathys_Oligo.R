
library("SingleCellExperiment")
library("tidyverse")
library("here")

sce <- readRDS(here("external-data", "AD_Knowledge_Portal", "Oligodendrocytes.rds"))
class(sce)

sce <- Seurat::as.SingleCellExperiment(sce)

colData(sce)

length((table(sce$projid))) 
#427
table(sce$cell_type_high_resolution)
# Oli 
# 645142

sce$barcode <- colnames(sce)
sce$barcode_end <- gsub("^[ACTG]+-","", sce$barcode)

pd <- as.data.frame(colData(sce))

pd |> count(projid, barcode_end) |> count(projid) |> count(n)
#   n  nn
# 1 1 378
# 2 2  48
# 3 3   1

#### dononr/clinical data ####
donor_info <- read_csv(here("external-data", "AD_Knowledge_Portal","ROSMAP_clinical.csv"),
                       col_types = cols(
                           projid = col_character()
                       ))

## fix leading 0 in projid
table(unique(sce$projid) %in% donor_info$projid)
# FALSE  TRUE 
# 13   414
missing_ids <- unique(sce$projid)[!unique(sce$projid) %in% donor_info$projid]

sce$projid <- as.character(sce$projid)
sce$projid[sce$projid %in% missing_ids] <- paste0("0", sce$projid[sce$projid %in% missing_ids])

all(unique(sce$projid) %in% donor_info$projid)

## filter donor info
donor_info_filter <- donor_info |>
    filter(projid %in% unique(sce$projid)) |>
    mutate(Age = ifelse(age_death == "90+", 90, as.numeric(age_death)))

summary(donor_info_filter$Age)
table(donor_info_filter$apoe_genotype)

table(donor_info_filter$braaksc)

table(donor_info_filter$braaksc, donor_info_filter$apoe_genotype)

table(donor_info_filter$ceradsc)
table(donor_info_filter$cogdx)
table(donor_info_filter$race)

all(!duplicated(donor_info_filter$individualID))

colData(sce)[sce$projid == "11409232",]

## example
donor_info_filter |> filter(projid == "11409232") |> select(individualID)

#### biospecimen data ####
# 
# biospec_data <- read_csv(here("external-data", "AD_Knowledge_Portal", "ROSMAP_biospecimen_metadata.csv")) |>
#     filter(individualID %in% donor_info_filter$individualID,
#            organ == "brain")
# 
# colData(sce)
# 
# biospec_data |> count(individualID, tissue)
# 
# table(donor_info_filter$individualID %in% biospec_data$individualID)


#### sn metadata ####

sn_metadata <- read_csv(here("external-data", "AD_Knowledge_Portal", "ROSMAP_assay_scrnaSeq_metadata.csv"))

table(sn_metadata$specimenID)
table(sn_metadata$libraryBatch)


#### multi-region data ####
## Mathys 2024, multiregion n=238 

## file syn52408579
mit_cell_labels <- read_delim(here("external-data", "AD_Knowledge_Portal","MIT_ROSMAP_Multiomics","all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv"),
                              col_types = cols(
                                  projid = col_character()
                              ))

all(unique(mit_cell_labels$projid) %in% donor_info$projid)

mit_cell_labels |> count(region, major.celltype)
mit_cell_labels |> count(projid, region)

# # A tibble: 283 × 3
# projid   region     n
# <chr>    <chr>  <int>
#     1 10291856 AG      4594
# 2 10291856 EC      3189
# 3 10291856 HC      7212
# 4 10291856 MT      5369

mit_cell_labels |> filter(major.celltype == "Oli") |> count(region)
# region     n
# <chr>  <int>
# 1 AG     57256
# 2 EC     69214
# 3 HC     78633
# 4 MT     48396
# 5 PFC    63192
# 6 TH     92265


#### hdf5 data ####
library("zellkonverter")

sce <- readH5AD(here("external-data", "AD_Knowledge_Portal","MIT_ROSMAP_Multiomics", "all_brain_regions_filt_preprocessed_scanpy_fullmatrix.h5ad"))

#### synapser package ? ####
## package install fails
# library(synapser) 
# synLogin(authToken="YOUR_TOKEN_HERE") 
# 
# # Obtain a pointer and download the data 
# syn52408579 <- synGet(entity='syn52408579') 

