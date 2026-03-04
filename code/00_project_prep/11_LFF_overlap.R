## Louise Huuki-Myers, Feb 2026
## compile all LFF project brains 

library("tidyverse")
library("here")

data_dir <- here("processed-data", "00_project_prep", "11_LFF_overlap")
if(!dir.exists(data_dir)) dir.create(data_dir)

select_cols <- c("BrNum", "Age", "Sex", "Ancestry", "APOE", "Rin")

## load data
erc_sample_info <- read.csv(here('processed-data', '00_project_prep', '05_pathology', 'sample_taupathy.csv')) |>
    select(all_of(select_cols)) |>
    mutate(ERC = TRUE,
           ERC_qc = BrNum != "Br1289")

sample_notes <- tibble(BrNum = "Br1289", note = "Wrong sample info")

## LC
lc_BrNum <- read_csv("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/00_project_prep/01_get_online_metadata/metadata_visium_plan.csv") |>
    transmute(Visium_slide = paste0(Visium_Slide, "_", Visium_subslide), 
             BrNum)

lc_sample_info <- read_csv("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/02_build_spe/sample_info.csv") |>
    left_join(lc_BrNum) |>
    select(all_of(select_cols)) |>
    mutate(LC = TRUE,
           APOE = gsub(" ", "", APOE),
           LC_qc = !BrNum %in% c("Br5517","Br5276","Br5712")) |>
    unique()

sample_notes <- sample_notes |> add_row(BrNum = c("Br5517","Br5276","Br5712"), note = "No LC in dissection")

## NbM
NbM_sample_info <- read_csv("/dcs05/lieber/lcolladotor/RAL_nbm_LIBD4265/RAL_nbm_APOE/processed-data/00_project_prep/RAL_nbm_APOE-sample_info.csv") |>
    select(all_of(select_cols)) |>
    mutate(NbM = TRUE)


LFF_sample_info <- erc_sample_info |>
    full_join(lc_sample_info, by = join_by(BrNum, Age, Sex, Ancestry, APOE, Rin)) |>
    full_join(NbM_sample_info, by = join_by(BrNum, Age, Sex, Ancestry, APOE, Rin)) |>
    replace_na(list(ERC = FALSE, LC = FALSE, NbM = FALSE)) |>
    left_join(sample_notes)

write_csv(LFF_sample_info, file = here(data_dir, "LFF_erc_lc_nbm_sample_info.csv"))

LFF_sample_info |>
    count(ERC_qc, LC_qc, NbM) |>
    arrange()

# ERC_qc LC_qc   NbM  n
# 1    TRUE  TRUE FALSE 11
# 2    TRUE  TRUE  TRUE 10
# 3      NA    NA  TRUE  9
# 4      NA  TRUE FALSE  7
# 5    TRUE    NA  TRUE  4
# 6    TRUE FALSE  TRUE  2
# 7    TRUE    NA FALSE  2
# 8      NA  TRUE  TRUE  2
# 9   FALSE    NA FALSE  1
# 10   TRUE FALSE FALSE  1
