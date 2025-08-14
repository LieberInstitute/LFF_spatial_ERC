## Louise Huuki-Myers, Aug 2025
## extract and format Grubman DE data

#### Set Up ####

library("tidyverse")
library("here")
library("sessioninfo")
library("readxl")

#### DE data ####

# sheets <- excel_sheets(here("external-data", "Grubman2019", "Grubman2019_SuppTables.xlsx"))
# 
# sheets_st4 <- sheets[grepl("Supplementary Table 4", sheets)]
# sheets_st4 <- sheets_st4[sheets_st4 != "Supplementary Table 4f"]
# 
# Grubman_DE_data <- map_dfr(sheets_st4,  ~read_excel(here("external-data", "Grubman2019", "Grubman2019_SuppTables.xlsx"), sheet = .x, skip =6)) |>
#     select(-`...1`) 
# 
# Grubman_DE_data |> count(cellType)
# Grubman_DE_data |> count(cellType)

# DE data downloaded from http://adsn.ddnetbio.com "AD vs Control cells within each cell type". LFC cutoff 0.5 (minimum option) FDR cutoff 0.05.

grubman_fn <- list.files(here("external-data", "Grubman2019"), pattern = ".csv", full.names = TRUE)
names(grubman_fn) <- gsub("Grubman2019_DEG_|.csv","",basename(grubman_fn))

Grubman_DE_data <- map2_dfr(grubman_fn, names(grubman_fn), ~read_csv(.x, show_col_types = FALSE) |> mutate(cell_type = .y))

# Grubman_DE_data |> count(cell_type)
# cell_type     n 
# <chr>     <int>.  # Extended data Fig. 3 values
# 1 Astro       705 #700
# 2 Endo        379
# 3 Micro       245 #193
# 4 Neuron      960 #285  - inhib 285 - excit 285
# 5 OPC         451 #320
# 6 Oligo       218 #228

# Grubman_DE_data |> arrange(-FDR)
# Grubman_DE_data |> arrange(abs(LFC))
