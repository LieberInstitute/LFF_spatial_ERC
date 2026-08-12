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
colors_path = here(
    'processed-data', '00_project_prep', 'cell_type_colors.V2.Rdata'
)
plot_dir = here('plots', '24_xenium_liana', '10_mediation_pairs')
n_samples_sig = 4

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

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

top_diff_density = function(
        filtered_df, full_df, colors, plot_path, n = 5, label_LR = TRUE
    ) {
    top_pairs = filtered_df |>
        slice_max(abs_lr_mean_diff, n = n, with_ties = FALSE)

    if (label_LR) {
        top_pairs = top_pairs |>
            mutate(
                label = sprintf(
                    '%s->%s\n%s->%s',
                    source, target, ligand_complex, receptor_complex
                )
            )
    } else {
        top_pairs = top_pairs |>
            mutate(label = sprintf('%s->%s', source, target))
    }

    pair_colors = colors[top_pairs$source]

    dens = density(full_df$lr_mean_diff, na.rm = TRUE)
    y_max = max(dens$y)

    # Assign y positions linearly spaced top-to-bottom, ordered by descending x
    y_positions = seq(0.9 * y_max, 0.2 * y_max, length.out = nrow(top_pairs))
    top_pairs = top_pairs |>
        arrange(desc(lr_mean_diff)) |>
        mutate(y_pos = y_positions)

    p = ggplot(full_df, aes(x = lr_mean_diff)) +
        geom_density(fill = 'steelblue', alpha = 0.5) +
        geom_vline(
            data = top_pairs,
            aes(xintercept = lr_mean_diff, color = source),
            linetype = 'dashed', linewidth = 1
        ) +
        ggrepel::geom_label_repel(
            data = top_pairs,
            aes(x = lr_mean_diff, y = y_pos, label = label, color = source),
            size = 6, inherit.aes = FALSE,
            direction = 'x', segment.size = 0.4, min.segment.length = 0
        ) +
        coord_cartesian(ylim = c(0, y_max)) +
        scale_color_manual(values = pair_colors) +
        labs(x = 'LR Mean % Difference', y = 'Density') +
        theme_bw(base_size = 18) +
        theme(legend.position = 'none')

    pdf(plot_path)
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
load(colors_path)

################################################################################
#   Astros as receiver; FZD8 as receptor
################################################################################

#-------------------------------------------------------------------------------
#   Most-different-pair distribution plot
#-------------------------------------------------------------------------------

score_comp_df |>
    filter(grepl('^Astro', target), receptor_complex == 'FZD8') |>
    top_diff_density(
        score_comp_df, cell_type_colors$anno,
        file.path(plot_dir, 'astro_FZD8_difference_distribution.pdf'),
        label_LR = FALSE
    )

#-------------------------------------------------------------------------------
#   Global scores across all samples
#-------------------------------------------------------------------------------

p = score_global_df |>
    #   The only ligand in any cell-type pair is IGFBP4
    filter(grepl('^Astro', target), receptor_complex == 'FZD8') |>
    slice_max(lr_mean, n = 5, with_ties = FALSE) |>
    arrange(desc(lr_mean)) |>
    mutate(
        cell_types = paste(source, '->', target),
        cell_types = factor(cell_types, levels = cell_types)
    ) |>
    ggplot(aes(x = cell_types, y = lr_mean, fill = source)) +
        geom_col() +
        theme_bw(base_size = 20) +
        guides(fill = 'none') +
        labs(x = 'Source -> Target', y = 'Mean LR Score') +
        scale_fill_manual(values = cell_type_colors$anno) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pdf(file.path(plot_dir, 'astro_FZD8_bar.pdf'), width = 5)
print(p)
dev.off()

score_global_df |>
    filter(target == 'Astro.2', receptor_complex == 'FZD8') |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    pair_density(
        score_global_df, file.path(plot_dir, 'astro2_FZD8_distribution.pdf')
    )

#-------------------------------------------------------------------------------
#   E4+ vs E2+ comparison
#-------------------------------------------------------------------------------

plot_df = score_comp_df |>
    filter(grepl('^Astro', target), receptor_complex == 'FZD8') |>
    mutate(cell_types = paste(source, '->', target))

p = ggplot(plot_df, aes(x = `lr_mean_E2+`, y = `lr_mean_E4+`)) +
    geom_point(color = 'grey70') +
    geom_point(
        data = plot_df |>
            slice_max(abs(lr_mean_diff), n = 5, with_ties = FALSE),
        aes(color = source),
        size = 3
    ) +
    geom_text_repel(
        data = plot_df |>
            slice_max(abs(lr_mean_diff), n = 5, with_ties = FALSE),
        aes(label = cell_types, color = source),
        size = 6
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
    theme_bw(base_size = 18) +
    guides(color = 'none') +
    labs(x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)') +
    scale_color_manual(values = cell_type_colors$anno) +
    scale_x_log10() +
    scale_y_log10()
pdf(file.path(plot_dir, 'astro_FZD8_scatter_difference.pdf'))
print(p)
dev.off()

message(
    sprintf(
        'Astro FZD8 plot only has one pair: %s->FZD8',
        unique(plot_df$ligand_complex)
    )
)

################################################################################
#   Oligo.3 as receiver; ERBB3 as receptor
################################################################################

#-------------------------------------------------------------------------------
#   Most-different-pair distribution plot
#-------------------------------------------------------------------------------

score_comp_df |>
    filter(target == 'Oligo.3', receptor_complex == 'ERBB3') |>
    top_diff_density(
        score_comp_df, cell_type_colors$anno,
        file.path(plot_dir, 'oligo3_ERBB3_difference_distribution.pdf'),
        label_LR = TRUE
    )

#-------------------------------------------------------------------------------
#   Global scores across all samples
#-------------------------------------------------------------------------------

p = score_global_df |>
    #   Only two ligands for these criteria
    filter(target == 'Oligo.3', receptor_complex == 'ERBB3') |>
    mutate(lr_pair = paste0(ligand_complex, '->', receptor_complex)) |>
    group_by(lr_pair) |>
    slice_max(lr_mean, n = 5, with_ties = FALSE) |>
    arrange(desc(lr_mean), .by_group = TRUE) |>
    mutate(source_plot = paste(source, lr_pair, sep = '___')) |>
    ungroup() |>
    mutate(source_plot = factor(source_plot, levels = unique(source_plot))) |>
    ggplot(aes(x = source_plot, y = lr_mean, fill = source)) +
        geom_col() +
        facet_wrap(~lr_pair, scales = 'free_x') +
        theme_bw(base_size = 20) +
        guides(fill = 'none') +
        labs(x = 'Source Cell Type', y = 'Mean LR Score') +
        scale_fill_manual(values = cell_type_colors$anno) +
        scale_x_discrete(labels = function(x) sub('___.*', '', x)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pdf(file.path(plot_dir, 'oligo3_ERBB3_bar.pdf'))
print(p)
dev.off()

score_global_df |>
    filter(target == 'Oligo.3', receptor_complex == 'ERBB3') |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    pair_density(
        score_global_df, file.path(plot_dir, 'oligo3_ERBB3_distribution.pdf')
    )

#-------------------------------------------------------------------------------
#   E4+ vs E2+ comparison
#-------------------------------------------------------------------------------

plot_df = score_comp_df |>
    filter(target == 'Oligo.3', receptor_complex == 'ERBB3') |>
    mutate(
        cell_types = paste(source, '->', target),
        lr_pair = paste0(ligand_complex, '->', receptor_complex)
    )

p = ggplot(
        plot_df, aes(x = `lr_mean_E2+`, y = `lr_mean_E4+`, shape = lr_pair)
    ) +
    geom_point(color = 'grey70') +
    geom_point(
        data = plot_df |>
            slice_max(abs(lr_mean_diff), n = 5, with_ties = FALSE),
        aes(color = source),
        size = 3
    ) +
    geom_text_repel(
        data = plot_df |>
            slice_max(abs(lr_mean_diff), n = 5, with_ties = FALSE),
        aes(label = cell_types, color = source),
        size = 6
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
    theme_bw(base_size = 20) +
    guides(color = 'none') +
    labs(
        x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)',
        shape = 'L->R'
    ) +
    scale_color_manual(values = cell_type_colors$anno) +
    scale_x_log10() +
    scale_y_log10()
pdf(file.path(plot_dir, 'oligo3_ERBB3_scatter_difference.pdf'), width = 10)
print(p)
dev.off()

################################################################################
#   APOE as ligand and any pair of cell types
################################################################################

#-------------------------------------------------------------------------------
#   Most-different-pair distribution plot
#-------------------------------------------------------------------------------

#   Only pair is APOE->TREM2
score_comp_df |>
    filter(ligand_complex == 'APOE') |>
    top_diff_density(
        score_comp_df, cell_type_colors$anno,
        file.path(plot_dir, 'APOE_difference_distribution.pdf'),
        label_LR = FALSE
    )

#-------------------------------------------------------------------------------
#   Global scores across all samples
#-------------------------------------------------------------------------------

p = score_global_df |>
    #   Only two ligands for these criteria
    filter(ligand_complex == 'APOE') |>
    mutate(lr_pair = paste0(ligand_complex, '->', receptor_complex)) |>
    group_by(lr_pair) |>
    slice_max(lr_mean, n = 5, with_ties = FALSE) |>
    arrange(desc(lr_mean), .by_group = TRUE) |>
    mutate(source_plot = paste(source, lr_pair, sep = '___')) |>
    ungroup() |>
    mutate(source_plot = factor(source_plot, levels = unique(source_plot))) |>
    ggplot(aes(x = source_plot, y = lr_mean, fill = source)) +
        geom_col() +
        facet_wrap(~lr_pair, scales = 'free_x') +
        theme_bw(base_size = 20) +
        guides(fill = 'none') +
        labs(x = 'Source Cell Type', y = 'Mean LR Score') +
        scale_fill_manual(values = cell_type_colors$anno) +
        scale_x_discrete(labels = function(x) sub('___.*', '', x)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pdf(file.path(plot_dir, 'APOE_bar.pdf'))
print(p)
dev.off()

score_global_df |>
    filter(ligand_complex == 'APOE') |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    pair_density(
        score_global_df, file.path(plot_dir, 'APOE_distribution.pdf')
    )

#-------------------------------------------------------------------------------
#   E4+ vs E2+ comparison
#-------------------------------------------------------------------------------

plot_df = score_comp_df |>
    filter(ligand_complex == 'APOE') |>
    mutate(
        cell_types = paste(source, '->', target),
        lr_pair = paste0(ligand_complex, '->', receptor_complex)
    )

p = ggplot(
        plot_df, aes(x = `lr_mean_E2+`, y = `lr_mean_E4+`, shape = lr_pair)
    ) +
    geom_point(color = 'grey70') +
    geom_point(
        data = plot_df |>
            slice_max(abs(lr_mean_diff), n = 5, with_ties = FALSE),
        aes(color = source),
        size = 3
    ) +
    geom_text_repel(
        data = plot_df |>
            slice_max(abs(lr_mean_diff), n = 5, with_ties = FALSE),
        aes(label = cell_types, color = source),
        size = 6
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
    theme_bw(base_size = 20) +
    guides(color = 'none') +
    labs(
        x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)',
        shape = 'L->R'
    ) +
    scale_color_manual(values = cell_type_colors$anno) +
    scale_x_log10() +
    scale_y_log10()
pdf(file.path(plot_dir, 'APOE_scatter_difference.pdf'), width = 10)
print(p)
dev.off()

session_info()
