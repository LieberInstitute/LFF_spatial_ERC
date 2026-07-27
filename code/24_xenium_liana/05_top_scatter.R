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
plot_dir = here('plots', '24_xenium_liana', '05_top_scatter')
n_samples_sig = 3

dir.create(plot_dir, showWarnings = FALSE)

score_df = read_csv(score_path, show_col_types = FALSE)

#   APOE is involved in a bunch of interactions, but not between Oligo and Astro
message('Cell types involved in APOE interactions:')
score_df |>
    filter(ligand_complex == 'APOE') |>
    mutate(
        source = sub('\\..*', '', source), target = sub('\\..*', '', target)
    ) |>
    distinct(source, target, APOE_carrier) |>
    print(n = Inf)

oligo_df = score_df |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    group_by(APOE_carrier, ligand_complex, receptor_complex) |>
    slice_max(lr_mean, with_ties = FALSE, n = 1) |>
    group_by(APOE_carrier) |>
    slice_max(lr_mean, with_ties = FALSE, n = 10) |>
    mutate(
        lr_pair = paste0(ligand_complex, '->', receptor_complex),
        cell_types = paste(source, '->', target)
    )

p = ggplot(oligo_df, aes(x = lr_mean, y = lr_specificity, color = cell_types)) +
    geom_point() +
    facet_wrap(~APOE_carrier) +
    geom_text_repel(aes(label = lr_pair), size = 5) +
    theme_bw(base_size = 20) +
    labs(x = 'Mean LR Score', y = 'LR Specificity', color = 'Source Cell Type')
pdf(file.path(plot_dir, 'oligo3_top_scatter_individual.pdf'), width = 12, height = 5)
print(p)
dev.off()

p = read_csv(unfiltered_path, show_col_types = FALSE) |>
    group_by(source, target, ligand_complex, receptor_complex, APOE_carrier) |>
    summarize(
        is_sig = sum(pval < 0.05) >= n_samples_sig,
        lr_mean = mean(lr_mean), # not necessarily among significant samples
    ) |>
    ungroup() |>
    pivot_wider(names_from = APOE_carrier, values_from = c(is_sig, lr_mean)) |>
    mutate(
        lr_mean_diff = 200 * (`lr_mean_E4+` - `lr_mean_E2+`) / (`lr_mean_E4+` + `lr_mean_E2+`),
        higher_in = ifelse(lr_mean_diff > 0, 'Higher in E4+', 'Higher in E2+'),
        lr_mean_diff = abs(lr_mean_diff)
    ) |>
    filter(ifelse(higher_in == 'Higher in E4+', `is_sig_E4+`, `is_sig_E2+`)) |>
    group_by(higher_in, ligand_complex, receptor_complex) |>
    mutate(lr_specificity = lr_mean_diff / max(lr_mean_diff)) |>
    ungroup() |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    group_by(higher_in, ligand_complex, receptor_complex) |>
    slice_max(lr_mean_diff, with_ties = FALSE, n = 1) |>
    group_by(higher_in) |>
    slice_max(lr_mean_diff, with_ties = FALSE, n = 10) |>
    mutate(
        lr_pair = paste0(ligand_complex, '->', receptor_complex),
        cell_types = paste(source, '->', target)
    ) |>
    ggplot(aes(x = lr_mean_diff, y = lr_specificity, color = cell_types)) +
        geom_point() +
        facet_wrap(~higher_in) +
        geom_text_repel(aes(label = lr_pair), size = 5) +
        theme_bw(base_size = 20) +
        labs(
            x = '% Diff in LR Score', y = 'LR Specificity',
            color = 'Source Cell Type'
        )
pdf(file.path(plot_dir, 'oligo3_top_scatter_difference.pdf'), width = 12, height = 5)
print(p)
dev.off()

p = read_csv(unfiltered_path, show_col_types = FALSE) |>
    group_by(source, target, ligand_complex, receptor_complex, APOE_carrier) |>
    summarize(lr_mean = mean(lr_mean)) |>
    ungroup() |>
    filter((source == 'Oligo.3') | (target == 'Oligo.3'), source != target) |>
    pivot_wider(names_from = APOE_carrier, values_from = lr_mean) |>
    ggplot(aes(x = `E2+`, y = `E4+`)) +
        geom_point(alpha = 0.3) +
        geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
        geom_smooth(method = 'lm', color = 'blue', se = FALSE) +
        theme_bw(base_size = 15) +
        labs(x = 'Mean LR Score (E2+)', y = 'Mean LR Score (E4+)') +
        scale_x_log10() +
        scale_y_log10()
pdf(file.path(plot_dir, 'oligo3_all_pairs.pdf'))
print(p)
dev.off()

spon1_app_value = oligo_df |>
    filter(
        source == 'Astro.5', ligand_complex == 'SPON1',
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
