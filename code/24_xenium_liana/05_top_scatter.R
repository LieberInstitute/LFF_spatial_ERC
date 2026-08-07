library(here)
library(tidyverse)
library(ggrepel)
library(sessioninfo)

score_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_summary.csv'
)
unfiltered_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_unfiltered.csv'
)
out_path = here(
    'processed-data', '24_xenium_liana', '05_top_scatter',
    'pair_level_comparison.csv'
)
plot_dir = here('plots', '24_xenium_liana', '05_top_scatter')
colors_path = here(
    'processed-data', '00_project_prep', 'cell_type_colors.V2.Rdata'
)
n_samples_sig = 3

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

individual_scatter = function(score_df, plot_path) {
    top_df = score_df |>
        group_by(APOE_carrier) |>
        slice_max(lr_mean, with_ties = FALSE, n = 10) |>
        mutate(
            lr_pair = paste0(ligand_complex, '->', receptor_complex),
            cell_types = paste(source, '->', target)
        )

    p = ggplot(top_df, aes(x = lr_mean, y = lr_specificity, color = cell_types)) +
        geom_point() +
        facet_wrap(~APOE_carrier) +
        geom_text_repel(aes(label = lr_pair), size = 5) +
        theme_bw(base_size = 20) +
        labs(x = 'Mean LR Score', y = 'LR Specificity', color = 'Source Cell Type')
    pdf(plot_path, width = 12, height = 5)
    print(p)
    dev.off()
}

comparison_scatter_manuscript = function(
    score_df, plot_path, color_var = cell_types, color_label = 'Other Cell Type',
    color_values = NULL
) {
    top_df = score_df |>
        group_by(higher_in) |>
        slice_max(lr_mean_diff, with_ties = FALSE, n = 10) |>
        mutate(
            lr_pair = paste0(ligand_complex, '->', receptor_complex),
            cell_types = paste(source, '->', target),
            other_cell_type = ifelse(source == 'Oligo.3', target, source),
            facet = ifelse(
                source == 'Oligo.3', 'Oligo.3 as Source', 'Oligo.3 as Target'
            )
        )

    p = ggplot(
        top_df,
        aes(x = lr_mean_diff, y = lr_specificity, color = {{ color_var }})
    ) +
        geom_point() +
        facet_grid(facet ~ higher_in) +
        geom_text_repel(aes(label = lr_pair), size = 7) +
        theme_bw(base_size = 25) +
        labs(
            x = '% Diff in LR Score', y = 'LR Specificity', color = color_label
        )

    if (!is.null(color_values)) {
        p = p + scale_color_manual(values = color_values)
    }

    pdf(plot_path, width = 12, height = 8)
    print(p)
    dev.off()
}

################################################################################
#   Load data
################################################################################

score_df = read_csv(score_path, show_col_types = FALSE)
unfiltered_sig_df = read_csv(unfiltered_path, show_col_types = FALSE) |>
    mutate(log_lr_mean_pseudo = log10(lr_mean + 1e-6)) |>
    group_by(source, target, ligand_complex, receptor_complex) |>
    summarize(
        n_e2 = sum(APOE_carrier == 'E2+'),
        n_e4 = sum(APOE_carrier == 'E4+'),
        n_nonzero_e2 = sum(APOE_carrier == 'E2+' & lr_mean > 0),
        n_nonzero_e4 = sum(APOE_carrier == 'E4+' & lr_mean > 0),
        `lr_mean_E2+` = mean(lr_mean[APOE_carrier == 'E2+'], na.rm = TRUE),
        `lr_mean_E4+` = mean(lr_mean[APOE_carrier == 'E4+'], na.rm = TRUE),
        p_value = {
            e2_vals = log_lr_mean_pseudo[APOE_carrier == 'E2+']
            e4_vals = log_lr_mean_pseudo[APOE_carrier == 'E4+']
            e2_var = var(e2_vals)
            e4_var = var(e4_vals)

            if (
                n_e2 < 2 || n_e4 < 2 ||
                    n_nonzero_e2 < 3 || n_nonzero_e4 < 3 ||
                    (isTRUE(e2_var == 0) && isTRUE(e4_var == 0))
            ) {
                NA_real_
            } else {
                t.test(e2_vals, e4_vals, alternative = 'two.sided')$p.value
            }
        },
        .groups = 'drop'
    ) |>
    filter(n_e2 >= 2, n_e4 >= 2, n_nonzero_e2 >= 3, n_nonzero_e4 >= 3) |>
    mutate(
        is_nom_sig = p_value < 0.01,
        avg_lr_mean = (`lr_mean_E2+` + `lr_mean_E4+`) / 2,
        lr_mean_diff = ifelse(
            avg_lr_mean == 0,
            0,
            200 * (`lr_mean_E4+` - `lr_mean_E2+`) /
                (`lr_mean_E4+` + `lr_mean_E2+`)
        ),
        higher_in = ifelse(lr_mean_diff > 0, 'Higher in E4+', 'Higher in E2+'),
        abs_lr_mean_diff = abs(lr_mean_diff)
    ) |>
    select(
        source, target, ligand_complex, receptor_complex, `lr_mean_E2+`,
        `lr_mean_E4+`, is_nom_sig, avg_lr_mean, lr_mean_diff, higher_in,
        abs_lr_mean_diff
    )

write_csv(unfiltered_sig_df, out_path)

load(colors_path)

#   APOE is involved in a bunch of interactions, but not between Oligo and Astro
message('Cell types involved in APOE interactions:')
score_df |>
    filter(ligand_complex == 'APOE') |>
    mutate(
        source = sub('\\..*', '', source), target = sub('\\..*', '', target)
    ) |>
    distinct(source, target, APOE_carrier) |>
    print(n = Inf)

################################################################################
#   Top interactions involving Oligo.3, by APOE carrier
################################################################################

#-------------------------------------------------------------------------------
#   Top interactions in E2+ and top in E4+
#-------------------------------------------------------------------------------

#   Oligo.3 and anything
score_df |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    group_by(APOE_carrier, ligand_complex, receptor_complex) |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    ungroup() |>
    individual_scatter(file.path(plot_dir, 'oligo3_top_scatter_individual.pdf'))

#   Astro -> Oligo.3
score_df |>
    filter(grepl('^Astro', source), target == 'Oligo.3') |>
    individual_scatter(
        file.path(plot_dir, 'oligo3_astro_top_scatter_individual.pdf')
    )

#-------------------------------------------------------------------------------
#   Top interactions higher in E2+ and higher in E4+
#-------------------------------------------------------------------------------

processed_df = unfiltered_sig_df |>
    filter(is_nom_sig) |>
    mutate(lr_mean_diff = abs_lr_mean_diff) |>
    group_by(higher_in, ligand_complex, receptor_complex) |>
    mutate(lr_specificity = lr_mean_diff / max(lr_mean_diff)) |>
    ungroup()

#   Oligo.3 and anything
processed_df |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    group_by(higher_in, ligand_complex, receptor_complex) |>
    slice_max(lr_mean_diff, with_ties = FALSE, n = 1) |>
    ungroup() |>
    comparison_scatter_manuscript(
        file.path(plot_dir, 'oligo3_top_scatter_difference.pdf'),
        color_var = other_cell_type,
        color_label = 'Other Cell Type',
        color_values = cell_type_colors$anno
    )

################################################################################
#   All pairs involving Oligo.3, comparing E2+ and E4+ (scatter plot)
################################################################################

p = unfiltered_sig_df |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    mutate(
        facet = ifelse(
            source == 'Oligo.3', 'Oligo.3 as Source', 'Oligo.3 as Target'
        ),
        other_cell_type = ifelse(source == 'Oligo.3', target, source),
        pair_label = ifelse(
            is_nom_sig,
            paste0(
                other_cell_type, ': ', ligand_complex, ' -> ', receptor_complex
            ),
            NA_character_
        )
    ) |>
    ggplot(
            aes(
                x = `lr_mean_E2+`, y = `lr_mean_E4+`, color = is_nom_sig,
                alpha = is_nom_sig
            )
        ) +
        geom_point() +
        geom_text_repel(aes(label = pair_label), size = 4, na.rm = TRUE) +
        geom_abline(
            slope = 1, intercept = 0, linetype = 'dashed', color = 'black'
        ) +
        facet_wrap(~facet) +
        theme_bw(base_size = 20) +
        scale_alpha_manual(values = c('FALSE' = 0.5, 'TRUE' = 1)) +
        guides(
            color = guide_legend(override.aes = list(alpha = 1)),
            alpha = 'none'
        ) +
        labs(
            x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)',
            color = 'p < 0.01'
        ) +
        scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = 'red3')) +
        scale_x_log10() +
        scale_y_log10()
pdf(file.path(plot_dir, 'oligo3_all_pairs.pdf'), width = 10, height = 5)
print(p)
dev.off()

################################################################################
#   All pairs comparing absolute and relative differences by carrier
################################################################################

p = unfiltered_sig_df |>
    ggplot(
            aes(
                x = avg_lr_mean, y = lr_mean_diff, color = is_nom_sig,
                alpha = is_nom_sig
            )
        ) +
        geom_point() +
        theme_bw(base_size = 18) +
        scale_color_manual(values = c('FALSE' = 'grey60', 'TRUE' = 'red3')) +
        scale_alpha_manual(values = c('FALSE' = 0.3, 'TRUE' = 1)) +
        guides(alpha = 'none') +
        labs(
            x = 'Mean LR Score (Average)', y = 'LR Score Difference (%)',
            color = 'p < 0.01'
        ) +
        scale_x_log10()
pdf(
    file.path(plot_dir, 'relative_and_absolute_diffs_all_pairs.pdf'),
    width = 8, height = 6
)
print(p)
dev.off()

################################################################################
#   Some ad-hoc requested plots
################################################################################

spon1_app_value = score_df |>
    filter(
        source == 'Astro.5', target == 'Oligo.3', ligand_complex == 'SPON1',
        receptor_complex == 'APP'
    ) |>
    pull(lr_mean)

p = ggplot(score_df, aes(x = lr_mean)) +
    geom_density(fill = 'steelblue', alpha = 0.5) +
    geom_vline(
        xintercept = spon1_app_value, linetype = 'dashed', color = 'red',
        linewidth = 1
    ) +
    annotate(
        'text', x = spon1_app_value, y = Inf,
        label = 'SPON1->APP\n(Astro.5->Oligo.3)', hjust = 1, vjust = 1.1,
        color = 'red', size = 6
    ) +
    labs(
        x = 'Mean LR Score', y = 'Density',
        title = 'Distribution of LR Scores\n(Significant Interactions)'
    ) +
    theme_bw(base_size = 20) +
    scale_x_log10()

pdf(file.path(plot_dir, 'lr_mean_distribution.pdf'), width = 9)
print(p)
dev.off()

p = score_df |>
    #   Basically Inhib.Vip <-> Astro or Oligo.3
    filter(
        grepl('^(Oligo.3|Astro|Inhib\\.Vip$)', source),
        grepl('^(Oligo.3|Astro|Inhib\\.Vip$)', target),
        source != target,
        (source == 'Inhib.Vip' | target == 'Inhib.Vip')
    ) |>
    #   We get a lot of results. Show only the source-target pair with the
    #   highest score for an LR pair
    group_by(ligand_complex, receptor_complex) |>
    slice_max(lr_mean, with_ties = FALSE) |>
    ungroup() |>
    #   Then show the top-10-scoring interactions
    arrange(desc(lr_mean)) |>
    slice_head(n = 10) |>
    mutate(
        lr_pair = paste0(ligand_complex, '->', receptor_complex),
        ct_pair = paste(source, '->', target)
    ) |>
    ggplot(aes(x = lr_mean, y = lr_specificity, color = ct_pair)) +
        geom_point() +
        geom_text_repel(aes(label = lr_pair), size = 4) +
        theme_bw(base_size = 20) +
        labs(
            x = 'Mean LR Score', y = 'LR Specificity',
            color = 'Cell Types'
        )
pdf(file.path(plot_dir, 'Inhib_Vip_top_scatter.pdf'), width = 9)
print(p)
dev.off()

session_info()
