library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(SpatialExperiment)

inflow_path = here(
    'processed-data', '24_xenium_liana', '08_prep_inflow', 'inflow_scores.csv'
)
spe_path = here(
    'processed-data', '21_Xenium', '13_xenium_bansky_embedding',
    'spe_xenium_bansky.qs2'
)
plot_dir = here('plots', '24_xenium_liana', '09_SpD_violin')

dir.create(plot_dir, showWarnings = FALSE)

spe = qs_read(spe_path)

inflow_df = read_csv(inflow_path, show_col_types = FALSE)
stopifnot(all(inflow_df$cell_id %in% colnames(spe)))

p = tibble(cell_id = colnames(spe), SpX = spe$SpX) |>
    inner_join(inflow_df, by = 'cell_id') |>
    group_by(SpX) |>
    filter(inflow_score > quantile(inflow_score, 0.99)) |>
    ungroup() |>
    ggplot(aes(x = SpX, y = inflow_score)) +
        geom_violin(fill = 'steelblue') +
        expand_limits(y = 0) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = 'Spatial Domain', y = 'Inflow Score')
pdf(file.path(plot_dir, 'SpD_violin.pdf'))
print(p)
dev.off()

session_info()
