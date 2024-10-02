## October 2024, Louise Huuki-Myers
## Check Pathology vs. genotype, create tile plot

library("tidyverse")
library("here")
library("sessioninfo")

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

plot_dir <- here("plots", "00_project_prep", "05_pathology")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## sample info
sample_info <- read_csv(here("processed-data", "02_build_spe",  "sample_info.csv")) |>
    mutate(APOE_carrier = ifelse(grepl("E2", APOE), "E2+", "E4+"))

## Read in data from RNAScope
taupathy <- read_csv(here("processed-data", "00_project_prep", "05_pathology", "Tau identification_KRM - Final Path Score.csv")) |>
    transmute(BrNum = gsub(" ", "", `Brain number`),
           taupathy = `Taupathy?` == "Yes")

taupathy |> count(taupathy)
# taupathy     n
# <lgl>    <int>
# 1 FALSE       24
# 2 TRUE         9

sample_info <- sample_info |>left_join(taupathy)

with(sample_info, table(APOE, taupathy))
#         taupathy
# APOE    FALSE TRUE
# E2/E2     4    2
# E2/E3     6    3
# E3/E4     9    2
# E4/E4     5    2

with(sample_info, table(Ancestry, taupathy))
#         taupathy
# Ancestry FALSE TRUE
# AA    13    4
# EA    11    5

with(sample_info, table(Sex, taupathy))
#     taupathy
# Sex FALSE TRUE
# F     7    2
# M    17    7

sample_info |>
    group_by(APOE_carrier) |>
    summarise(prop_tau = sum(taupathy)/n())
# APOE_carrier prop_tau
# <chr>           <dbl>
# 1 E2+             0.333
# 2 E4+             0.222

sample_info |>
    group_by(APOE) |>
    summarise(prop_tau = sum(taupathy)/n())
# APOE  prop_tau
# <chr>    <dbl>
# 1 E2/E2    0.333
# 2 E2/E3    0.333
# 3 E3/E4    0.182
# 4 E4/E4    0.286

## save data
sample_info |>
    select(BrNum, APOE, APOE_carrier, Ancestry, Sex, Age, Diagnosis, Anc_Afr, Anc_Eur, Rin, taupathy) |>
    write_csv(here("processed-data", "00_project_prep", "05_pathology","sample_taupathy.csv"))
