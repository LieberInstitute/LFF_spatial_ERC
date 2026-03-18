## Louise Huuki-Myers, March 2025
## Check Oligo marker genes from external data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ComplexHeatmap")

data_dir <- here("processed-data", "19_other_Oligo", "00_other_Oligo_markers")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "00_other_Oligo_markers")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Compile external marker gene data ####

# Jakel et al. (2019)
sheets <- readxl::excel_sheets(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"))

Oligo_stats_jakel <- map_dfr(sheets[2:10], ~readxl::read_xlsx(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"), sheet = .x) |>
                                 mutate(Jakel_Oligo = paste0("Jakel_", .x))) |>
    select(gene, 
           logFC = avg_logFC,
           Oligo = Jakel_Oligo,
           ) |>
    mutate(Dataset = "Jakel2019")

# Mathys et al. (2019)
## TODO

# Sadick et al. (2022)
Oligo_stats_sadick_stab3 <- readxl::read_xlsx(here("external-data", "Sadick2022","Sadick2022_STab3.xlsx"), sheet = "LEN_so_oligo_DEGs") |>
    mutate(Oligo = paste0("Sadick_S", LEN_so_oligo_cluster),
           Dataset = "Sadick2022") |>
    select(gene,
           logFC = avg_log2FC,
           Oligo,
           Dataset)


# Grubman et al. (2019)
grubman_subtypes_fn <- list.files(here("external-data", "Grubman2019"), pattern = "Oligo_", full.names = TRUE)
names(grubman_subtypes_fn) <- gsub("Grubman_Oligo_|.csv", "", basename(grubman_subtypes_fn))

grubman_subtypes_data <- map2_dfr(grubman_subtypes_fn, names(grubman_subtypes_fn), ~read.csv(.x) |> 
                                      mutate(g_cell_type = .y)) |>
    rename(gene = geneName)

# Siletti et al. (2023) Allen brain atlas
read_csv(here("external-data", "Siletti2023", "science.add7046_table_s5.csv")) |>
    rename(gene = `...1`, Siletti_stat = stat)