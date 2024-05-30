
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("rtracklayer")
library("lobstr")
library("tidyverse")
library("sessioninfo")

## Create output directories
dir_rdata <- here::here("processed-data", "02.1_spe_compare_trim")
if(!dir.exists(dir_rdata)) dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

#### Read in Sample Info ####
## check datatype, use factors when possible
sample_info <- tibble(sample_id = list.files(here("processed-data", "01_spaceranger")),
                      sample_path = here(list.files(here("processed-data", "01_spaceranger"), full.names = TRUE), "outs")) |>
    separate(sample_id, into = c("sample", "trimmed"), sep = "_u", remove = FALSE) |>
    mutate(trimmed = is.na(trimmed)) |>
    group_by(sample) |>
    mutate(any_trimmed = any(!trimmed),
           n = n())

## select only sample with trimmed and untrimmed samples
sample_info <- sample_info |> filter(any_trimmed) |> select(-any_trimmed, - n)

dim(sample_info)
# [1] 30  4

table(sample_info$trimmed)

## all files exist
stopifnot(all(file.exists(sample_info$sample_path)))
message("Processing data for ", sum(file.exists(sample_info$sample_path)), " samples...")

sample_info$sample_path[!file.exists(sample_info$sample_path)]

## updated to untrimmed
# list.files(here("processed-data", "01_spaceranger"))[!paste0(list.files(here("processed-data", "01_spaceranger"), full.names = TRUE),"/outs") %in% sample_info$sample_path]
# [1] "V13B23-363_A1" "V13B23-363_B1" "V13B23-363_C1" "V13B23-363_D1" "V13B23-364_A1" "V13B23-364_B1" "V13B23-364_C1"
# [8] "V13B23-364_D1" "V13B23-365_A1" "V13B23-365_B1" "V13B23-365_C1" "V13B23-365_D1" "V13B23-366_B1" "V13B23-366_C1"
# [15] "V13B23-366_D1"

## write csv for easy access
# write.csv(sample_info, here(dir_rdata, "sample_info.csv"), row.names = FALSE)


## Build basic SPE
message(Sys.time(), "- Starting read10xVisiumWrapper")
spe <- read10xVisiumWrapper(
  samples = sample_info$sample_path,
  sample_id = sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE
)
message(Sys.time(), "- Done read10xVisiumWrapper")

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

table(spe$sample_id)
table(spe$trimmed)

## use key as colnames to prevent duplicate barcodes
colnames(spe) <- spe$key

## Save with and without dropping spots outside of the tissue
spe_trim <- spe

## Size in Gb
message("Size of spe_trim:")
lobstr::obj_size(spe_trim)
# 5.56 GB

## save
message(Sys.time(), "- Saving")
saveRDS(spe_trim, here(dir_rdata, "spe_trim.rds"))

# slurmjobs::job_single('01_build_spe_trim', create_shell = TRUE, memory = '25G', command = "Rscript 01_build_spe_trim.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
