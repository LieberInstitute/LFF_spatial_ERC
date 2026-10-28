library(here)
library(tidyverse)
library(ggrepel)
library(sessioninfo)

comparison_path = here(
    'processed-data', '24_xenium_liana', '05_top_scatter',
    'pair_level_comparison.csv'
)
deg_path = here(
    'processed-data', '17_cell_cell_comm', '10_LR_DEG_check',
    'Liana_LR_DEGs_sn_fine.csv'
)
colors_path = here(
    'processed-data', '00_project_prep', 'cell_type_colors.V2.Rdata'
)
plot_dir = here('plots', '24_xenium_liana', '11_sn_DEG_overlap_plots')
pair_columns = c('source', 'target', 'ligand_complex', 'receptor_complex')
n_top = 5

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
load(colors_path)

################################################################################
#   Load data and restrict to sn DEG-overlapping directed interactions
################################################################################

deg_pairs = read_csv(deg_path, show_col_types = FALSE) |>
    distinct(across(all_of(pair_columns)))
stopifnot(!anyNA(deg_pairs))

#   Preserve the source/target context in which the ligand or receptor is a DEG.
unfiltered_sig_df = read_csv(comparison_path, show_col_types = FALSE) |>
    semi_join(deg_pairs, by = pair_columns)
stopifnot(nrow(unfiltered_sig_df) > 0)

################################################################################
#   Highlight the five largest absolute percent differences per Oligo.3 facet
################################################################################

plot_df = unfiltered_sig_df |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    mutate(
        facet = ifelse(
            source == 'Oligo.3', 'Oligo.3 as Source', 'Oligo.3 as Target'
        ),
        other_cell_type = ifelse(source == 'Oligo.3', target, source),
        pair_label = paste0(ligand_complex, ' -> ', receptor_complex)
    )

#   No overlapping pairs are nominally significant (p < 0.01) in these data,
#   so show the largest absolute percent differences instead.
top_pairs = plot_df |>
    group_by(facet) |>
    slice_max(abs(lr_mean_diff), n = n_top, with_ties = FALSE) |>
    ungroup()

p = ggplot(plot_df, aes(x = `lr_mean_E2+`, y = `lr_mean_E4+`)) +
    geom_point(color = 'grey70', size = 2) +
    geom_point(
        data = top_pairs,
        aes(color = other_cell_type),
        size = 3.5
    ) +
    geom_text_repel(
        data = top_pairs,
        aes(label = pair_label, color = other_cell_type),
        size = 5, max.overlaps = Inf, seed = 7381,
        box.padding = 0.6, force = 10, max.time = 5, show.legend = FALSE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
    facet_wrap(~facet) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'bottom') +
    guides(color = guide_legend(nrow = 2)) +
    labs(
        x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)',
        color = 'Other Cell Type'
    ) +
    scale_color_manual(values = cell_type_colors$anno) +
    scale_x_log10() +
    scale_y_log10()

pdf(file.path(plot_dir, 'oligo3_all_pairs.pdf'), width = 10, height = 6)
print(p)
dev.off()

session_info()
