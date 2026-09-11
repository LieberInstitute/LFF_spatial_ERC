#   Fix cell-type names in the main data sources

library(tidyverse)
library(here)
library(sessioninfo)
library(readxl)

cell_map_path = here(
    'processed-data', '04_snRNA-seq', '27.5_sn_neuron_check',
    'ERC_neuron_subcluster_details_annotated.xlsx'
)
summary_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_summary.csv'
)
unfiltered_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_unfiltered.csv'
)
out_dir = here(
    'processed-data', '24_xenium_liana', '03_2_update_cell_types'
)

dir.create(out_dir, showWarnings = FALSE)

cell_map_df = read_excel(cell_map_path) |>
    select(cell_type_anno, cell_type_update)
cell_map = setNames(cell_map_df$cell_type_update, cell_map_df$cell_type_anno)

summary_df = read_csv(summary_path, show_col_types = FALSE) |>
    mutate(
        source = ifelse(is.na(cell_map[source]), source, cell_map[source]),
        target = ifelse(is.na(cell_map[target]), target, cell_map[target])
    )

unfiltered_df = read_csv(unfiltered_path, show_col_types = FALSE) |>
    mutate(
        source = ifelse(is.na(cell_map[source]), source, cell_map[source]),
        target = ifelse(is.na(cell_map[target]), target, cell_map[target])
    )

write_csv(summary_df, file.path(out_dir, 'global_interactions_summary.csv'))
write_csv(
    unfiltered_df, file.path(out_dir, 'global_interactions_unfiltered.csv')
)

session_info()
