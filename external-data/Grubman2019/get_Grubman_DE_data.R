## Louise Huuki-Myers, Aug 2025
## extract and format Grubman DE data

#### Set Up ####

library("tidyverse")
library("here")
library("sessioninfo")
library("readxl")

#### DE data ####

sheets <- excel_sheets(here("external-data", "Grubman2019", "Grubman2019_SuppTables.xlsx"))

sheets_st4 <- sheets[grepl("Supplementary Table 4", sheets)]
sheets_st4 <- sheets_st4[sheets_st4 != "Supplementary Table 4f"]

Grubman_DE_data <- map_dfr(sheets_st4,  ~read_excel(here("external-data", "Grubman2019", "Grubman2019_SuppTables.xlsx"), sheet = .x, skip =6)) |>
    select(-`...1`) 

Grubman_DE_data |> count(cellType)


# Mathys_DE_data |> dplyr::count(cell_type, DE_type)
