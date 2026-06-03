library(here)
library(tidyverse)
library(ggrepel)
library(sessioninfo)

score_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_summary.csv'
)
plot_dir = here('plots', '24_xenium_liana', '05_top_scatter')

dir.create(plot_dir, showWarnings = FALSE)

score_df = read_csv(score_path, show_col_types = FALSE)

#   APOE is involved in a bunch of interactions, but not between Oligo and Astro
message('Cell types involved in APOE interactions:')
score_df |>
    filter(ligand_complex == 'APOE') |>
    mutate(
        source = sub('\\..*', '', source), target = sub('\\..*', '', target)
    ) |>
    distinct(source, target) |>
    print(n = Inf)

astro_oligo = score_df |>
    filter(grepl('^Astro', source), target == 'Oligo.3') |>
    mutate(lr_pair = paste0(ligand_complex, '->', receptor_complex))

p = ggplot(astro_oligo, aes(x = lr_mean, y = lr_specificity, color = source)) +
    geom_point() +
    geom_text_repel(aes(label = lr_pair), size = 4) +
    theme_bw(base_size = 20) +
    labs(x = 'Mean LR Score', y = 'LR Specificity', color = 'Source Cell Type')
pdf(file.path(plot_dir, 'astro_oligo3_top_scatter.pdf'), width = 9)
print(p)
dev.off()

spon1_app_value = astro_oligo |>
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

session_info()
