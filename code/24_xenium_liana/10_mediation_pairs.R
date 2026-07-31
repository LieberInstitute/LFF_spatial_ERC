library(here)
library(tidyverse)
library(sessioninfo)
library(ggrepel)

score_global_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_unfiltered.csv'
)
score_comp_path = here(
    'processed-data', '24_xenium_liana', '05_top_scatter',
    'pair_level_comparison.csv'
)
plot_dir = here('plots', '24_xenium_liana', '10_mediation_pairs')
n_samples_sig = 4

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

pair_scatter = function(score_df, plot_path) {
    p = score_df |>
        arrange(desc(lr_mean)) |>
        slice_head(n = 10) |>
        mutate(
            lr_pair = paste0(ligand_complex, '->', receptor_complex),
            cell_types = paste(source, '->', target)
        ) |>
        ggplot(aes(x = lr_mean, y = lr_specificity, color = cell_types)) +
            geom_point() +
            geom_text_repel(aes(label = lr_pair), size = 5) +
            theme_bw(base_size = 20) +
            labs(x = 'Mean LR Score', y = 'LR Specificity', color = 'Cell Types')
    pdf(plot_path, width = 9)
    print(p)
    dev.off()
}

pair_density = function(pair_df, score_df, plot_path) {
    pair_df = pair_df |>
        mutate(
            label = sprintf(
                '%s->%s\n%s->%s',
                ligand_complex, receptor_complex, source, target
            )
        )

    p = ggplot(score_df, aes(x = lr_mean)) +
        geom_density(fill = 'steelblue', alpha = 0.5) +
        geom_vline(
            xintercept = pair_df$lr_mean, linetype = 'dashed', color = 'red',
            linewidth = 1
        ) +
        annotate(
            'text', x = pair_df$lr_mean, y = Inf, label = pair_df$label, hjust = 1,
            vjust = 1.1, color = 'red', size = 6
        ) +
        labs(
            x = 'Mean LR Score', y = 'Density',
            title = 'Distribution of LR Scores\n(Significant Interactions)'
        ) +
        theme_bw(base_size = 20) +
        scale_x_log10()

    pdf(plot_path, width = 9)
    print(p)
    dev.off()
}

################################################################################
#   Process significant pairs
################################################################################

score_global_df = read_csv(score_global_path, show_col_types = FALSE) |>
    filter(pval < 0.05) |>
    group_by(source, target, ligand_complex, receptor_complex) |>
    filter(length(unique(sample_id)) >= n_samples_sig) |>
    summarize(lr_mean = mean(lr_mean)) |>
    group_by(ligand_complex, receptor_complex) |>
    mutate(lr_specificity = lr_mean / max(lr_mean)) |>
    ungroup()

score_comp_df = read_csv(score_comp_path, show_col_types = FALSE)

################################################################################
#   Astros as receiver; FZD8 as receptor
################################################################################

#-------------------------------------------------------------------------------
#   Global scores across all samples
#-------------------------------------------------------------------------------

score_global_df |>
    #   The only ligand in any cell-type pair is IGFBP4
    filter(grepl('^Astro', target), receptor_complex == 'FZD8') |>
    pair_scatter(file.path(plot_dir, 'astro_FZD8_scatter.pdf'))

score_global_df |>
    filter(target == 'Astro.2', receptor_complex == 'FZD8') |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    pair_density(
        score_global_df, file.path(plot_dir, 'astro2_FZD8_distribution.pdf')
    )

#-------------------------------------------------------------------------------
#   E4+ vs E2+ comparison
#-------------------------------------------------------------------------------

score_comp_df |>
    filter(grepl('^Astro', target), receptor_complex == 'FZD8') |>
    ggplot(aes(x = `lr_mean_E2+`, y = `lr_mean_E4+`)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
        theme_bw(base_size = 20) +
        labs(
            x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)'
        ) +
        scale_x_log10() +
        scale_y_log10()
    pair_scatter(file.path(plot_dir, 'astro_FZD8_scatter_difference.pdf'))
################################################################################
#   Oligo.3 as receiver; ERBB3 as receptor
################################################################################

score_df |>
    #   Only two ligands for these criteria
    filter('Oligo.3' == target, receptor_complex == 'ERBB3') |>
    pair_scatter(file.path(plot_dir, 'oligo3_ERBB3_scatter.pdf'))

score_df |>
    filter(target == 'Oligo.3', receptor_complex == 'ERBB3') |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    pair_density(score_df, file.path(plot_dir, 'oligo3_ERBB3_distribution.pdf'))

session_info()
