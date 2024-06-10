
# cd /dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/

# remotes::install_github("drighelli/SpatialExperiment")
# remotes::install_github("LieberInstitute/spatialLIBD")

library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
library("dplyr")

## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
if(!dir.exists(dir_rdata)) dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "02_build_spe", "01_build_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Read in Sample Info ####
message(Sys.time(), "- Read in Sample Info")
## check datatype, use factors when possible
sample_info <- read.csv(here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_visium_plan.csv")) |>
    filter(is.na(lc_note)) |>
    mutate(Visium_slide = paste0(`Visium.Slide..`, "_" ,`Visium_subslide`),
           APOE = gsub(", ", "/", APOE,),
           sample_path = here("processed-data", "01_spaceranger", paste0(Visium_slide,"_untrimmed"), ## if untrimmed file exist, select it
                              "outs"),
           Ancestry = gsub("CAUC", "EA", Ancestry)) |>
    select(sample_id = BrNum, Visium_slide, APOE, Ancestry, Sex, Age, Diagnosis, Rin, sample_path) |>
    mutate(sample_path = ifelse(file.exists(sample_path),
                                sample_path, 
                                gsub("_untrimmed", "", sample_path)),
           base_path = gsub("^.*?/(V.*?)/outs","\\1", sample_path),
           BrNum = sample_id)

## all files exist
stopifnot(all(file.exists(sample_info$sample_path)))
message("Processing data for ", sum(file.exists(sample_info$sample_path)), " samples...")

# sample_info$sample_path[!file.exists(sample_info$sample_path)]


## write csv for easy access
write.csv(sample_info, here(dir_rdata, "sample_info.csv"), row.names = FALSE)

## Define some info for the samples

## Define the donor info using information from - moving away from providing data by hand, did double check
## https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/raw-data/sample_info/Visium_HPC_Round1_20220113_Master_ADR.xlsx
## https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/raw-data/sample_info/Visium_HPC_Round2_20220223_Master_ADR.xlsx
# donor_info <- data.frame(
# sample_id = c("V13Y24-343_A1","V13Y24-343_B1","V13Y24-343_C1","V13Y24-343_D1","V13Y24-344_A1","V13Y24-344_B1","V13Y24-344_C1","V13Y24-344_D1","V13Y24-342_A1","V13Y24-342_B1","V13Y24-342_C1","V13Y24-342_D1","V13Y24-340_A1","V13Y24-340_B1","V13Y24-340_C1","V13Y24-340_D1","V13B23-363_A1","V13B23-363_B1","V13B23-363_C1","V13B23-363_D1","V13B23-364_A1","V13B23-364_B1","V13B23-364_C1","V13B23-364_D1","V13B23-365_A1","V13B23-365_B1","V13B23-365_C1","V13B23-365_D1","V13B23-366_B1","V13B23-366_C1","V13B23-366_D1"),
#   age = c(51.63, 51.45, 29.95, 59.86, 41.44, 60.83, 62.70, 46.53, 42.39, 48.75, 50.08, 61.34, 63.98, 54.88, 67.75, 31.31, 42.19, 55.88, 45.3, 51.11, 59.98, 43.67, 68.38, 58.19, 51.73, 50.2, 54.43, 48.69, 57.1, 60.84, 50.73),
#   sex = c("M","F","F","M", "M", "M", "M", "F", "F", "F", "M", "M", "M", "F", "M", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F","M", "M", "F"),
#   race = c("EA/CAUC", "AA", "AA", "EA/CAUC", "EA/CAUC", "AA", "AA", "AA", "AA", "AA", "EA/CAU", "AA", "EA/CAUC", "AA", "EA/CAUC", "EA/CAUC", "AA", "EA/CAUC", "AA", "EA/CAUC", "AA", "EA/CAUC", "AA", "AA", "EA/CAU", "AA", "EA/CAU", "AA", "EA/CAU", "EA/CAU", "AA"),
#   diagnosis = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control"),
#   rin = c(8.5, 8.6, 9.3, 7.4, 7.3, 8.7, 7.2, 8.5, 7.1, 9.2, 8.2, 8.5, 9.4, 8, 6.8, 6.6, 9, 8.5, 8.4, 8.6, 9.3, 5.2, 9, 8.7, 7.1, 8.3, 8.5, 5.2, 7.7, 7.4, 7), # fix rin information for V13B23-363_C1 sample
#   apoe = c("E2/E3", "E4/E4", "E2/E2", "E3/E4", "E2/E2", "E2/E3","E3/E4","E4/E4", "E3/E4", "E2/E3", "E4/E4", "E3/E4", "E3/E4", "E4/E4", "E2/E2", "E2/E3", "E4/E4", "E3/E4", "E2/E2", "E2/E3", "E3/E4", "E2/E2", "E2/E3", "E4/E4", "E2/E3", "E3/E4", "E4/E4", "E2/E2", "E3/E4", "E2/E3", "E3/E4")
# )
# 
# identical(donor_info$sample_id, sample_info$sample_id)
# identical(donor_info$age, sample_info$Age)
# identical(donor_info$sex, sample_info$Sex)
# identical(donor_info$race, sample_info$Ancestry) # "EA/CAUC" diff
# identical(donor_info$diagnosis, sample_info$Diagnosis)
# identical(donor_info$rin, sample_info$Rin)
# identical(donor_info$apoe, sample_info$APOE)

#### Build basic SPE ####
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

# [1] "2024-04-16 13:54:16 EDT"
# 2024-04-16 13:54:22.104989 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2024-04-16 14:07:31.501473 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2024-04-16 14:07:59.115353 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2024-04-16 14:08:18.078968 rtracklayer::import: reading the reference GTF file
# 2024-04-16 14:09:36.644289 adding gene information to the SPE object
# 2024-04-16 14:09:37.673669 adding information used by spatialLIBD
# [1] "2024-04-16 14:09:54 EDT"

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

## use key as colnames to prevent duplicate barcodes
colnames(spe) <- spe$key
any(duplicated(colnames(spe)))

## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

#### Read in cell counts and segmentation results ####
message(Sys.time(), " - Add cell segmentation")
sample_info |> 
    mutate(spot_count = file.exists(here(sample_path, "spatial","tissue_spot_counts.csv"))) |>
    select(sample_id, base_path, spot_count)

message("Existing spot count files: ", 
        sum(file.exists(here(sample_info$sample_path, "spatial","tissue_spot_counts.csv"))),
        "/",
        nrow(sample_info))

segmentations_list <-
    lapply(sample_info$sample_path, function(path) {
        file <-
            here(path,
                 "spatial",
                 "tissue_spot_counts.csv"
            )
        if (!file.exists(file)) {
            return(NULL)
        }
        x <- read.csv(file)
        x$Visium_slide <- gsub("_untrimmed","", gsub("^.*?/(V.*?)/outs","\\1", path))
        x$key <- paste0(x$barcode, "_", sample_info$sample_id[which(path == sample_info$sample_path)])
        return(x)
    })

## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <-
#     Reduce(function(...) {
#         merge(..., all = TRUE)
#     }, segmentations_list[lengths(segmentations_list) > 0])

segmentations <- do.call("rbind", segmentations_list)
dim(segmentations) 
# [1] 154752     11

# ## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
    segmentations[segmentation_match, -which(
        colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key", "Visium_slide")
    )]

head(segmentation_info) 
colData(spe) <- cbind(colData(spe), segmentation_info)

#### QC ####
message(Sys.time(), "- Add Quality Control metrics")
## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
message("Genes with no expression: ", length(no_expr), " (", round(length(no_expr) / nrow(spe) * 100, 3), "%)")
# [1]  6126
# [1] 16.73725
spe <- spe[-no_expr, ]

spe_raw <- add_qc_metrics(spe)

table(spe_raw$in_tissue, spe_raw$scran_discard)

table(spe_raw$sample_id, spe_raw$scran_low_lib_size_edge)

spe_qc_notes <- as.data.frame(colData(spe_raw)) |>
    filter(in_tissue) |>
    group_by(sample_id) |>
    summarize(n = n(),
              min_sum_umi = min(sum_umi),
              median_sum_umi = median(sum_umi),
              max_sum_umi = max(sum_umi),
              min_sum_gene = min(sum_gene),
              median_sum_gene = median(sum_gene),
              max_sum_gene = max(sum_gene),
              min_expr_chrM_ratio = min(expr_chrM_ratio),
              median_expr_chrM_ratio = median(expr_chrM_ratio),
              max_expr_chrM_ratio = max(expr_chrM_ratio),
              n_low_lib_edge = sum(scran_low_lib_size_edge == "TRUE", na.rm = TRUE),
              n_scran_discard = sum(scran_discard == "TRUE", na.rm = TRUE))

write.csv(spe_qc_notes, file = here(dir_rdata, "spe_qc_notes.csv"), row.names = FALSE)

## record cutoffs
# as.data.frame(colData(spe_raw)) |>
#     filter(in_tissue, as.logical(scran_low_lib_size)) |>
#     group_by(sample_id) |>
#     summarize(max_sum_umi = max(sum_umi), 
#               n = n()) |>
#     arrange(max_sum_umi) |>
#     print(n = 31)

## Size in Gb
message("Size of spe_raw:")
lobstr::obj_size(spe_raw)
# 5.59 GB

message(Sys.time(), " - Saving Data")
saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

message(Sys.time(), " - Prelim plots")
vis_grid_clus(
    spe = spe_raw,
    clustervar = "in_tissue",
    pdf = here::here(plot_dir, "spe_grid_in_tissue.pdf"),
    sort_clust = FALSE,
    point_size = 1
    # colors = c("in" = "grey90", "out" = "orange")
)

vis_grid_clus(
    spe = spe_raw,
    clustervar = "scran_low_lib_size_edge",
    pdf = here::here(plot_dir, "spe_grid_scran_low_lib_size_edge.pdf"),
    sort_clust = FALSE,
    point_size = 1
    # colors = c("in" = "grey90", "out" = "orange")
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
