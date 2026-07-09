# Create a tibble mapping future array tasks to subsetting behavior: each array
# task will potentially subset data by spatial domain, or subset to Oligo.3
# cells near astrocytes with high/low APOE expression

library(here)
library(tidyverse)
library(sessioninfo)
library(qs2)
library(SpatialExperiment)

spe_path = here(
    'processed-data', '21_Xenium', '13_xenium_bansky_embedding',
    'spe_xenium_bansky.qs2'
)
out_path = here(
    'processed-data', '14_MOFA', '07_build_tasks', 'task_map.csv'
)

dir.create(dirname(out_path), showWarnings = FALSE)

spe = qs_read(spe_path)

rbind(
        tibble(domain = levels(spe$SpX), astro_group = 'all'),
        tibble(domain = 'all', astro_group = c('high', 'low', 'all'))
    ) |>
    mutate(task_id = row_number()) |>
    write_csv(out_path)

session_info()
