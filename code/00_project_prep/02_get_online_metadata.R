

library("tidyverse")
library("googlesheets4")
library("sessioninfo")
library("here")

data_dir <- here("processed-data", "00_project_prep", "02_get_online_metadata")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## Study
metadata_study <- read_sheet("https://docs.google.com/spreadsheets/d/1mHEIBhN7kckInOipyi2ozqtYKLPeAbw1SL6KskWRCLM/edit#gid=0",
                       sheet = "Study",
                       skip = 1) |>
    rename(BrNum  = `Brain #`) |>
    filter(grepl("^Br", BrNum))

metadata_study |> count(`LC: Right/Left/ Both, slab`)

write_csv(metadata_study, file = here(data_dir, "metadata_study.csv"))

## Visium List
metadata_visium_list <- read_sheet("https://docs.google.com/spreadsheets/d/1mHEIBhN7kckInOipyi2ozqtYKLPeAbw1SL6KskWRCLM/edit#gid=0",
                             sheet = "Visium list") |>
    rename(BrNum  = `Brain #`, visium_note = `...10`) |>
    mutate(redisscet = grepl("re-dis", BrNum),
           Visium = !(grepl("No visium", visium_note) | grepl("Stopped|stopped|drop", Notes)),
           BrNum = gsub(" \\(re-dis\\)", "", BrNum)) 

metadata_visium_list |> count(Visium)
# Visium     n
# <lgl>  <int>
# 1 FALSE      5
# 2 TRUE      30
metadata_visium_list |> filter(Visium)|> count(APOE)
# APOE       n
# <chr>  <int>
# 1 E2, E2     5
# 2 E2, E3     8
# 3 E3, E4    10
# 4 E4, E4     7
any(duplicated(metadata_visium_list$BrNum))

write_csv(metadata_visium_list, file = here(data_dir, "metadata_visium_list.csv"))

## Visium Plan 
metadata_visium_plan <- read_sheet("https://docs.google.com/spreadsheets/d/1mHEIBhN7kckInOipyi2ozqtYKLPeAbw1SL6KskWRCLM/edit#gid=0",
                                   sheet = "Visium plan")|>
    filter(!is.na(`Visium slide 1`),
           Donors != "Donors") |>
    rename(BrNum  = `Donors`, 
           Visium_subslide  = `Visium slide 1`,
           lc_note = `...13`) |>
    mutate(redisscet = grepl("re-dis", BrNum),
           BrNum = gsub("\\(re-dis\\)", "", BrNum),
           remade_libaries = grepl("remade libraries", Date),
           Date = ifelse(remade_libaries, NA, Date),
           Rin = as.double(Rin),
           Age = as.double(Age),
           SBox = as.integer(SBox)
           )  |> 
    fill(Date, `Visium Slide #`) |>
    add_column(slide_index = rep(1:8, each = 4))

colnames(metadata_visium_plan)
# [1] "Date"            "Visium Slide #"  "Visium_subslide" "BrNum"           "Genotype"        "Age"             "Sex"            
# [8] "Ancestry"        "Diagnosis"       "Rin"             "APOE"            "SBox"            "lc_note"         "redisscet"      
# [15] "remade_libaries" "slide_index" 

write_csv(metadata_visium_plan, file = here(data_dir, "metadata_visium_plan.csv"))

## snRNA-seq plan
metadata_sn_plan <- read_sheet("https://docs.google.com/spreadsheets/d/1mHEIBhN7kckInOipyi2ozqtYKLPeAbw1SL6KskWRCLM/edit#gid=0",
                                   sheet = "snRNA-seq plan",
                               range = "A:O")|>
    filter(!is.na(Donors),
           Donors != "Donors") |>
    rename(BrNum  = `Donors`, 
           snRNA_complete  = `snRNA-seq round 1`) |>
    mutate(redisscet = grepl("re-dis", BrNum),
           BrNum = gsub("\\(re-dis\\)", "", BrNum),
           snRNA_complete = snRNA_complete == "Complete",
           APOE = gsub(",E", ", E", APOE)) |>
    replace_na(list(snRNA_complete = FALSE))

metadata_sn_plan |> count(snRNA_complete)
metadata_sn_plan |> count(APOE)

write_csv(metadata_sn_plan, file = here(data_dir, "metadata_sn_plan.csv"))

## RNAscope
metadata_rnascope <- read_sheet("https://docs.google.com/spreadsheets/d/1mHEIBhN7kckInOipyi2ozqtYKLPeAbw1SL6KskWRCLM/edit#gid=0",
                               sheet = "ERC_RNAscope_panels")|>
    filter(!is.na(Gene)) |>
    fill(Panel)

write_csv(metadata_rnascope, file = here(data_dir, "metadata_rnascope.csv"))

# slurmjobs::job_single('02_get_online_metadata', create_shell = TRUE, memory = '5G', command = "Rscript 02_get_online_metadata.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
