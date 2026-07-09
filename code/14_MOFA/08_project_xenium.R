library(here)
library(tidyverse)
library(sessioninfo)
library(qs2)
library(SpatialExperiment)

spe_path = here(
    'processed-data', '21_Xenium', '13_xenium_bansky_embedding',
    'spe_xenium_bansky.qs2'
)
astro_path = here(
    'processed-data', '21_Xenium', '20_xenium_Oligo3_Astro',
    'Oligo3_Astro_neighbor_df.Rdata'
)
task_path = here(
    'processed-data', '14_MOFA', '07_build_tasks', 'task_map.csv'
)

this_task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe = qs_read(spe_path)

astro_df = load(astro_path) |>
    get() |>
    dplyr::rename(cell_id = reference_barcode) |>
    select(cell_id, APOE_level)

spe$astro_group = tibble(cell_id = spe$cell_id) |>
    left_join(astro_df, by = 'cell_id') |>
    pull(APOE_level)

task_df = read_csv(task_path, show_col_types = FALSE) |>
    filter(task_id == this_task_id)

#   Filter Oligo.3 as applicable
if (task_df$astro_group != 'all') {
    spe = spe[
        ,
        (spe$astro_group == task_df$astro_group) &
        (spe$cell_type_anno != 'Oligo.3')
    ]
}

#   Filter spatial domain as applicable
if (task_df$domain != 'all') {
    spe = spe[, spe$SpX == task_df$domain]
}
