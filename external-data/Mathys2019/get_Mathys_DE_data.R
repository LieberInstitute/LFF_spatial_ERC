## Louise Huuki-Myers, Aug 2025
## extract and format Mathys DE data

#### Set Up ####

library("tidyverse")
library("here")
library("sessioninfo")
library("readxl")

#### DE data ####

# excel_sheets(here("external-data", "Mathys2019", "Mathys2019_SuppTable6.xlsx"))

cell_type_sheets <- c("Ex","In", "Ast","Oli","Opc","Mic")

data_range <- c(`no_v_path` = "A2:I20000", no_v_early = "L2:T20000", early_v_late = "W2:AE20000")

Mathys_DE_data <- map_dfr(cell_type_sheets, function(sheet){
    map2_dfr(data_range, names(data_range), 
             ~read_excel(here("external-data", "Mathys2019", "Mathys2019_SuppTable6.xlsx"), sheet, range = .x) |>
                 dplyr::rename(gene_name = `...1`) |>
                 dplyr::filter(!is.na(gene_name)) |>
                 dplyr::mutate(DE_type = .y,
                               cell_type = sheet,
                               .before = 1) |>
                 dplyr::mutate(MixedModel.z = as.double(MixedModel.z),
                               MixedModel.p = as.double(MixedModel.p))
             )
})

# Mathys_DE_data |> dplyr::count(cell_type, DE_type)
