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

clean_boxplot = function(plot_df, plot_path) {
    # Compute whisker bounds per domain to set y-axis limits
    whisker_bounds = plot_df |>
        group_by(SpX) |>
        summarise(
            lower = quantile(inflow_score, 0.25) - 1.5 * IQR(inflow_score),
            upper = quantile(inflow_score, 0.75) + 1.5 * IQR(inflow_score),
            .groups = 'drop'
        )

    y_max = max(whisker_bounds$upper)

    p = ggplot(plot_df, aes(x = SpX, y = inflow_score)) +
        geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = c(0, y_max)) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = 'Spatial Domain', y = 'Inflow Score')
    pdf(plot_path)
    print(p)
    dev.off()
}

spe = qs_read(spe_path)

inflow_df = read_csv(inflow_path, show_col_types = FALSE)
stopifnot(all(inflow_df$cell_id %in% colnames(spe)))

plot_df = tibble(
        cell_id = colnames(spe), SpX = spe$SpX, cell_type = spe$cell_type_anno
    ) |>
    inner_join(inflow_df, by = 'cell_id')

clean_boxplot(plot_df, file.path(plot_dir, 'SpD_boxplot_all_cell_types.pdf'))
clean_boxplot(
    plot_df |>
        filter(cell_type == 'Astro.3'),
    file.path(plot_dir, 'SpD_boxplot_astro_3.pdf')
)

session_info()
