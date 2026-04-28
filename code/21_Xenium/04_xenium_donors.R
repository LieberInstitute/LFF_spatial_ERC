## Louise Huuki-Myers, Feb 2026
## select 12 samples for Xenium experiment

#### Set up ####
library("tidyverse")
library("here")

data_dir <- here("processed-data", "21_Xenium", "04_xenium_donors")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "04_xenium_donors")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


## 21 common Brains LC - ERC
common_BrNum <- readLines(here("processed-data", "00_project_prep", "ERC_LC_common_BrNum.txt"))

## LC samples with fibroblast 
fibroblast <- c("Br6297","Br6538","Br6423","Br5529","Br5426","Br5415","Br5276","Br5368","Br6222","Br6263","Br5941","Br6098","Br6119")

## rep sections
rep_section <- read_csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "rep_section.csv")) |>
    select(BrNum = sample_id, rep_section)

## sample notes 

## Exclude
## Mislabled sample Br1289
## Outliers in DE LOO Br6263, Br5529, Br6423, Br3974

flag_samples <- tibble(BrNum = c("Br6263", "Br5529", "Br6423", "Br3974"), flag = "DE LOO") |>
    add_row(BrNum = "Br1691", flag = "MOFA")

sample_info <- sample_info <- read_csv(here("processed-data", "00_project_prep", "05_pathology","sample_taupathy.csv")) |>
    mutate(LC = BrNum %in% common_BrNum,
           fibroblast = BrNum %in% fibroblast) |>
    left_join(rep_section) |>
    filter(BrNum != "Br1289")  |>
    left_join(flag_samples)


sample_info |>
    count(fibroblast, rep_section, LC, is.na(flag))

sample_info |>
    count(LC, fibroblast)

sample_info |>
    count(is.na(flag), APOE_carrier, Ancestry)

sample_info |>
    filter(fibroblast) |>
    count(APOE_carrier, Ancestry)

sample_info |> 
    filter(rep_section) |> 
    group_by(APOE_carrier, Ancestry) |>
    count()

sample_info |>
    filter(LC) |> 
    group_by(APOE_carrier, Ancestry) |>
    count()


sample_info_score <- sample_info |> 
    filter(is.na(flag)) |>
    rowwise() |>
    mutate(score = sum(rep_section, LC, 2*fibroblast)) |>
    select(BrNum, APOE_carrier, APOE, Ancestry, Sex, Age, LC, rep_section, fibroblast, score) 

sample_info_score |>
    group_by(APOE_carrier, Ancestry) |> 
    count(score) |>
    pivot_wider(names_from = "score", values_from = "n")

# APOE_carrier Ancestry   `0`   `1`   `2`   `3`
# <chr>        <chr>    <int> <int> <int> <int>
# 1 E2+          AA           1     1     2*     1*
# 2 E2+          EA           3     2     1*    1*
# 3 E4+          AA          NA     3     4*    NA
# 4 E4+          EA           3     2*     1*    NA

sample_info_score |>
    arrange(-score) |>
    filter(APOE_carrier == "E4+", Ancestry == "EA")

sample_info |>
    filter(APOE_carrier == "E4+", Ancestry == "EA")
    
(sample_picks <- sample_info_score |>
        arrange(-score) |>
        group_by(APOE_carrier, Ancestry) |>
        mutate(rank = row_number()) |>
        filter(rank <= 4) |>
        arrange(APOE_carrier, Ancestry, rank))


#    BrNum  APOE_carrier Ancestry Sex     Age LC    rep_section fibroblast score  rank
#    <chr>  <chr>        <chr>    <chr> <dbl> <lgl> <lgl>       <lgl>      <dbl> <int>
#  1 Br5415 E2+          AA       M      60.8 TRUE  TRUE        TRUE           4     1
#  2 Br5426 E2+          AA       F      48.8 TRUE  FALSE       TRUE           3     2
#  3 Br1556 E2+          AA       F      30.0 TRUE  TRUE        FALSE          2     3 (1 M, 2F)
#  4 Br5854 E2+          AA       M      68.4 TRUE  FALSE       FALSE          1     4 * good
#  5 Br6538 E2+          EA       M      60.8 TRUE  TRUE        TRUE           4     1
#  6 Br2582 E2+          EA       M      41.4 TRUE  TRUE        FALSE          2     2
#  7 Br5367 E2+          EA       F      31.3 TRUE  FALSE       FALSE          1     3 (2 M, 1 F)
#  8 Br5634 E2+          EA       M      51.1 TRUE  FALSE       FALSE          1     4 * 3/4 tie
#  9 Br5941 E4+          AA       M      42.2 TRUE  FALSE       TRUE           3     1
# 10 Br6098 E4+          AA       M      50.2 TRUE  FALSE       TRUE           3     2
# 11 Br1039 E4+          AA       F      51.4 TRUE  TRUE        FALSE          2     3 (2 M, 1 F)
# 12 Br5460 E4+          AA       M      61.3 TRUE  TRUE        FALSE          2     4 * 3/4 tie
# 13 Br5276 E4+          EA       M      59.9 FALSE FALSE       TRUE           2     1
# 14 Br1706 E4+          EA       M      50.1 TRUE  TRUE        FALSE          2     2
# 15 Br5517 E4+          EA       M      64.0 FALSE TRUE        FALSE          1     3 (3 M)
# 16 Br5599 E4+          EA       M      55.9 FALSE FALSE       FALSE          0     4 *good

write_csv(sample_picks, file = here(data_dir, "xenium_donor_picks.csv"))

(sample_picks <- sample_info_score |>
        arrange(-score) |>
        filter(APOE_carrier == "E4+", Ancestry == "EA") |>
        mutate(rank = row_number()) |>
        arrange(APOE_carrier, Ancestry, rank))


sample_picks_confirmed <- readxl::read_xlsx(here(data_dir, "xenium_donor_picks_tissue_inventory.xlsx")) 

# sample_picks_confirmed |> dplyr::count(`Tissue left?`)

#### plot Visium data for selected samples ####

sample_info <- read_csv(here("processed-data", "21_Xenium", "08_xenium_QC_normalize", "ERC_Xenium_QC_summary.csv"))


sample_info |> select(BrNum) |> write_csv(here(data_dir, "Xenium_rotation.csv"))

library("spatialLIBD")

message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

vis_clus_plots <- vis_grid_clus(
    spe = spe[, spe$BrNum %in% sample_info$BrNum],
    clustervar = "SpD",
    colors = metadata(spe)$SpD_colors,
    sort_clust = FALSE,
    return_plots = FALSE,
    point_size = 2,
    guide_point_size = 3,
    pdf_file = here(plot_dir, "Visium_SpD_xenium_samples.pdf")
)




