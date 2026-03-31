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

#### Load external marker gene data ####
Oligo_stats <- list()

# Jakel et al. (2019)
sheets <- readxl::excel_sheets(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"))

Oligo_stats$Jakel <- map_dfr(sheets[2:10], ~readxl::read_xlsx(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"), sheet = .x) |>
                                 mutate(Jakel_Oligo = paste0("Jakel_", .x))) |>
    select(gene, 
           logFC = avg_logFC,
           Oligo = Jakel_Oligo,
           ) |>
    mutate(Dataset = "Jakel2019")

# Mathys et al. (2019)
## TODO

# Sadick et al. (2022)
Oligo_stats$Sadick <- readxl::read_xlsx(here("external-data", "Sadick2022","Sadick2022_STab3.xlsx"), sheet = "LEN_so_oligo_DEGs") |>
    mutate(Oligo = paste0("Sadick_S", LEN_so_oligo_cluster),
           Dataset = "Sadick2022") |>
    select(gene,
           logFC = avg_log2FC,
           Oligo,
           Dataset)


# Grubman et al. (2019)
grubman_subtypes_fn <- list.files(here("external-data", "Grubman2019"), pattern = "Oligo_", full.names = TRUE)
names(grubman_subtypes_fn) <- gsub("Grubman_Oligo_|.csv", "", basename(grubman_subtypes_fn))

Oligo_stats$Grubman <- map2_dfr(grubman_subtypes_fn, names(grubman_subtypes_fn), ~read.csv(.x) |> 
                                      mutate(Oligo = paste0("Grubman_", .y),
                                             Dataset = "Grubman2019")
                                ) |>
                                    select(gene = geneName,
                                           logFC = LFC,
                                           Oligo,
                                           Dataset)

Oligo_stats$Grubman |> count(Oligo)

# Siletti et al. (2023) Allen brain atlas
Oligo_stats$Siletti <- read_csv(here("external-data", "Siletti2023", "science.add7046_table_s5.csv")) |>
    rename(gene = `...1`, Siletti_stat = stat) |>
    mutate(Oligo = "Siletti_OPALIN",
           Dataset = "Siletti2023") |>
    select(gene,
           logFC = log2FoldChange,
           Oligo,
           Dataset)

#### Compile external marker gene data ####

Oligo_stats_all <- do.call("bind_rows", Oligo_stats)

Oligo_stats_all |> filter(gene %in% c("RBFOX1", "RASGRF1", "OPALIN")) |> arrange(gene)

Oligo_stats_all |> count(Dataset, Oligo) 
gene_count <- Oligo_stats_all |> group_by(gene) |> arrange(-n)
gene_count |> filter(n > 2)
# |>
#     dplyr::summarise(n = dplyr::n(), .by = c(Oligo)) |>
#     dplyr::filter(n > 1L) 

Oligo_stats_wide <- Oligo_stats_all |> 
    filter(gene %in% (gene_count |> filter(n>2) |> pull(gene)))|>
    select(-Dataset) |>
    pivot_wider(names_from = "Oligo", values_from = "logFC") |>
    column_to_rownames("gene") |>
    as.matrix()

Heatmap(Oligo_stats_wide, cluster_rows = FALSE, cluster_columns = FALSE)

Oligo_stats_wide[1:5, 1:5]

Oligo_stats_wide_lgl <- is.na(Oligo_stats_wide)
mode(Oligo_stats_wide_lgl) <- "numeric"
# Oligo_stats_wide_lgl[1:5, 1:5]
Heatmap(Oligo_stats_wide_lgl, cluster_rows = TRUE, cluster_columns = TRUE)
    
Oligo_stats_wide_na0 <- Oligo_stats_wide
Oligo_stats_wide_na0[is.na(Oligo_stats_wide_na0)] <- 0
Heatmap(Oligo_stats_wide_na0, cluster_rows = TRUE, cluster_columns = TRUE)

