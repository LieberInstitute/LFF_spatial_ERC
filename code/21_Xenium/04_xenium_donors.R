## Louise Huuki-Myers, Feb 2026
## select 12 samples for Xenium experiment

#### Set up ####
library("tidyverse")
library("here")


## 21 common Brains LC - ERC
common_BrNum <- readLines(here("processed-data", "00_project_prep", "ERC_LC_common_BrNum.txt"))

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
    mutate(LC = BrNum %in% common_BrNum) |>
    left_join(rep_section) |>
    filter(BrNum != "Br1289")  |>
    left_join(flag_samples)


sample_info |>
    count(rep_section, LC, is.na(flag))

sample_info |>
    count(is.na(flag), APOE_carrier, Ancestry)

sample_info |> 
    filter(rep_section) |> 
    group_by(APOE_carrier, Ancestry) |>
    count()

sample_info |>
    filter(LC) |> 
    group_by(APOE_carrier, Ancestry) |>
    count()


sample_info |> 
    filter(is.na(flag)) |>
    rowwise() |>
    mutate(score = sum(LC, rep_section)) |>
    select(BrNum, APOE_carrier, Ancestry, LC, rep_section, score) |>
    group_by(APOE_carrier, Ancestry) |> 
    count(score) |>
    pivot_wider(names_from = "score", values_from = "n")

#    APOE_carrier Ancestry   `0`   `1`   `2`
# <chr>        <chr>    <int> <int> <int>
# 1 E2+          AA           1     2     2
# 2 E2+          EA           3     2     2
# 3 E4+          AA          NA     5     2
# 4 E4+          EA           4     1     1
    
    filter(score > 0)
    slice(1:3)

    



