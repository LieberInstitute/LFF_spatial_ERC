library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggrepel)
library(sessioninfo)

score_global_path = here(
    'processed-data', '24_xenium_liana', '03_2_update_cell_types',
    'global_interactions_unfiltered.csv'
)
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
n_samples_sig = 3
pair_columns = c('source', 'target', 'ligand_complex', 'receptor_complex')

cell_type_levels = c(
    'Astro.1', 'Astro.2', 'Astro.3', 'Astro.4', 'Astro.5', 'Macro', 'Micro.1',
    'Micro.2', 'Micro.3', 'Micro.4', 'Micro.5', 'OPC.1', 'OPC.2', 'OPC.3',
    'OPC.4', 'OPC.5', 'Oligo.1', 'Oligo.2', 'Oligo.3', 'Oligo.4', 'Oligo.5',
    'Vasc.Endo', 'Vasc.PC', 'Vasc.VLMC', 'Excit.L2', 'Excit.L2_5.mix',
    'Excit.L5.RELN', 'Excit.L5.RORB', 'Excit.L2_5.RELN', 'Excit.L5_6.NP',
    'Excit.L6.CT', 'Excit.L6b', 'Inhib.HTR3A', 'Inhib.LAMP5_LHX6',
    'Inhib.PVALB', 'Inhib.VIP', 'Inhib.Chandelier', 'Inhib.SST'
)
carrier_groups = c('E2+', 'E4+')

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
load(colors_path)

################################################################################
#   Functions from 04_02_global_score_heatmap.R
################################################################################

make_diff_matrix = function(diff_df) {
    diff_df |>
        select(source, target, symmetric_pct_diff) |>
        pivot_wider(
            names_from = target,
            values_from = symmetric_pct_diff,
            values_fill = 0
        ) |>
        column_to_rownames('source') |>
        as.matrix() |>
        _[cell_type_levels, cell_type_levels, drop = FALSE]
}

save_difference_heatmap = function(diff_df, out_path) {
    diff_mat = make_diff_matrix(diff_df)
    max_abs_diff = max(abs(diff_mat), na.rm = TRUE)
    max_abs_diff = if_else(max_abs_diff > 0, max_abs_diff, 1)
    annotation_colors = cell_type_colors$anno[cell_type_levels]
    cell_type_labels = factor(cell_type_levels, levels = cell_type_levels)

    ht = Heatmap(
        diff_mat,
        name = sprintf(
            'Sym. %% diff\n%s vs %s', carrier_groups[[2]], carrier_groups[[1]]
        ),
        col = colorRamp2(
            c(-max_abs_diff, 0, max_abs_diff),
            c('#3b4cc0', 'white', '#b40426')
        ),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_dend_side = 'right',
        row_names_side = 'left',
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 90,
        top_annotation = HeatmapAnnotation(
            Target = cell_type_labels,
            col = list(Target = annotation_colors),
            annotation_name_gp = gpar(fontsize = 9),
            show_legend = FALSE
        ),
        left_annotation = rowAnnotation(
            Source = cell_type_labels,
            col = list(Source = annotation_colors),
            annotation_name_gp = gpar(fontsize = 9),
            show_legend = FALSE
        ),
        heatmap_legend_param = list(title = 'Sym. % diff')
    )

    pdf(out_path, width = 10, height = 8)
    draw(
        ht,
        column_title = sprintf(
            'Symmetric percent difference in summed mean lr_mean (%s vs %s)',
            carrier_groups[[2]], carrier_groups[[1]]
        ),
        column_title_gp = gpar(fontsize = 14, fontface = 'bold')
    )
    dev.off()
}

################################################################################
#   Load data and restrict to sn DEG-overlapping directed interactions
################################################################################

deg_pairs = read_csv(deg_path, show_col_types = FALSE) |>
    distinct(across(all_of(pair_columns)))
stopifnot(!anyNA(deg_pairs))

#   Preserve the source/target context in which the ligand or receptor is a DEG.
score_global_df = read_csv(score_global_path, show_col_types = FALSE) |>
    semi_join(deg_pairs, by = pair_columns)

#   Reuse step 05's carrier means and nominal significance without retesting.
unfiltered_sig_df = read_csv(comparison_path, show_col_types = FALSE) |>
    semi_join(deg_pairs, by = pair_columns)
stopifnot(nrow(score_global_df) > 0, nrow(unfiltered_sig_df) > 0)

message(
    nrow(distinct(score_global_df, across(all_of(pair_columns)))), ' of ',
    nrow(deg_pairs), ' sn DEG interactions present in Xenium scores; ',
    nrow(unfiltered_sig_df), ' retained in step 05 comparisons.'
)

################################################################################
#   Difference heatmap, as in 04_02_global_score_heatmap.R
################################################################################

avg_lr = score_global_df |>
    filter(pval < 0.05) |>
    group_by(source, ligand_complex, receptor_complex, target, APOE_carrier) |>
    filter(n_distinct(sample_id) >= n_samples_sig) |>
    summarize(lr_mean = mean(lr_mean), .groups = 'drop')

heatmap_df = expand_grid(
        source = cell_type_levels,
        target = cell_type_levels,
        APOE_carrier = carrier_groups
    ) |>
    left_join(
        avg_lr |>
            group_by(source, target, APOE_carrier) |>
            summarize(lr_mean = sum(lr_mean), .groups = 'drop'),
        by = c('source', 'target', 'APOE_carrier')
    ) |>
    mutate(
        lr_mean = replace_na(lr_mean, 0),
        source = factor(source, levels = cell_type_levels),
        target = factor(target, levels = cell_type_levels),
        APOE_carrier = factor(APOE_carrier, levels = carrier_groups)
    )

diff_df = heatmap_df |>
    select(source, target, APOE_carrier, lr_mean) |>
    pivot_wider(
        names_from = APOE_carrier,
        values_from = lr_mean,
        values_fill = 0
    ) |>
    mutate(
        denominator = .data[[carrier_groups[[2]]]] + .data[[carrier_groups[[1]]]],
        symmetric_pct_diff = if_else(
            denominator == 0,
            0,
            200 * (.data[[carrier_groups[[2]]]] -
                .data[[carrier_groups[[1]]]]) / denominator
        )
    )

save_difference_heatmap(
    diff_df,
    file.path(
        plot_dir,
        sprintf(
            'source_target_heatmap_sig%s_APOE_difference_complexheatmap.pdf',
            n_samples_sig
        )
    )
)

################################################################################
#   All pairs involving Oligo.3, as in 05_top_scatter.R
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

session_info()
