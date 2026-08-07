library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(sessioninfo)

score_global_path = here(
    'processed-data', '24_xenium_liana', '04_global_score_heatmap',
    'global_interactions_unfiltered.csv'
)
colors_path = here(
    'processed-data', '00_project_prep', 'cell_type_colors.V2.Rdata'
)
plot_dir = here('plots', '24_xenium_liana', '04_02_global_score_heatmap')
n_samples_sig = 3

cell_type_levels = c(
    'Astro.1', 'Astro.2', 'Astro.3', 'Astro.4', 'Astro.5', 'Macro', 'Micro.1',
    'Micro.2', 'Micro.3', 'Micro.4', 'Micro.5', 'OPC.1', 'OPC.2', 'OPC.3',
    'OPC.4', 'OPC.5', 'Oligo.1', 'Oligo.2', 'Oligo.3', 'Oligo.4', 'Oligo.5',
    'Vasc.Endo', 'Vasc.PC', 'Vasc.VLMC', 'Excit.L2', 'Excit.L2_5.1',
    'Excit.L2_5.2', 'Excit.L5.1', 'Excit.L5.2', 'Excit.L5_6_NP', 'Excit.L6_CT',
    'Excit.L6b', 'Inhib.Pax6', 'Inhib.Lamp5_Lhx6', 'Inhib.Pvalb', 'Inhib.Vip',
    'Inhib.Chandelier', 'Inhib.Sst'
)
carrier_groups = c('E2+', 'E4+')
broad_cell_type_levels = c(
    'Astro', 'Macro', 'Micro', 'OPC', 'Oligo', 'Vasc', 'Excit', 'Inhib'
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
load(colors_path)

################################################################################
#   Functions
################################################################################

get_broad_cell_type = function(cell_type) {
    case_when(
        grepl('^Astro', cell_type) ~ 'Astro',
        grepl('^Macro', cell_type) ~ 'Macro',
        grepl('^Micro', cell_type) ~ 'Micro',
        grepl('^OPC', cell_type) ~ 'OPC',
        grepl('^Oligo', cell_type) ~ 'Oligo',
        grepl('^Vasc', cell_type) ~ 'Vasc',
        grepl('^Excit', cell_type) ~ 'Excit',
        grepl('^Inhib', cell_type) ~ 'Inhib',
        .default = NA_character_
    )
}

make_lr_matrix = function(df, carrier_group) {
    df |>
        filter(APOE_carrier == carrier_group) |>
        select(source, target, lr_mean) |>
        pivot_wider(
            names_from = target,
            values_from = lr_mean,
            values_fill = 0
        ) |>
        arrange(factor(source, levels = cell_type_levels)) |>
        column_to_rownames('source') |>
        as.matrix() |>
        _[, cell_type_levels, drop = FALSE]
}

make_source_target_heatmap = function(
    mat, name, title, fill_fun, show_heatmap_legend = TRUE
) {
    cell_type_broad = factor(
        get_broad_cell_type(cell_type_levels),
        levels = broad_cell_type_levels
    )
    cell_type_labels = factor(cell_type_levels, levels = cell_type_levels)
    annotation_colors = cell_type_colors$anno[cell_type_levels]

    Heatmap(
        mat,
        name = paste(name, title, sep = '_'),
        col = fill_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = cell_type_broad,
        column_split = cell_type_broad,
        row_order = cell_type_levels,
        column_order = cell_type_levels,
        row_gap = unit(2, 'mm'),
        column_gap = unit(2, 'mm'),
        row_title = NULL,
        column_title = title,
        column_title_gp = gpar(fontsize = 13, fontface = 'bold'),
        row_names_side = 'left',
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 90,
        show_heatmap_legend = show_heatmap_legend,
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
        heatmap_legend_param = list(title = name)
    )
}

save_faceted_global_heatmap = function(heatmap_df, out_path) {
    max_lr_mean = max(heatmap_df$lr_mean, na.rm = TRUE)
    fill_fun = colorRamp2(
        c(0, max_lr_mean * c(1e-6, 0.33, 0.66, 1)),
        c('lightgray', '#440154', '#31688e', '#35b779', '#fde725')
    )

    heatmaps = lapply(
        seq_along(carrier_groups),
        \(i) {
            carrier_group = carrier_groups[[i]]
            make_source_target_heatmap(
                make_lr_matrix(heatmap_df, carrier_group),
                name = 'Sum mean\nlr_mean',
                title = carrier_group,
                fill_fun = fill_fun,
                show_heatmap_legend = i == 1
            )
        }
    )

    pdf(out_path, width = 16, height = 8)
    draw(
        heatmaps[[1]] + heatmaps[[2]],
        column_title = sprintf(
            'Sum of mean lr_mean for interactions significant in >=%s samples',
            n_samples_sig
        ),
        column_title_gp = gpar(fontsize = 15, fontface = 'bold'),
        merge_legends = TRUE
    )
    dev.off()
}

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
        name = sprintf('Sym. %% diff\n%s vs %s', carrier_groups[[2]], carrier_groups[[1]]),
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
#   Global communication heatmaps
################################################################################

score_global_df = read_csv(score_global_path, show_col_types = FALSE)

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

save_faceted_global_heatmap(
    heatmap_df,
    file.path(
        plot_dir,
        sprintf('source_target_heatmap_sig%s_by_APOE_carrier_complexheatmap.pdf', n_samples_sig)
    )
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
            200 * (.data[[carrier_groups[[2]]]] - .data[[carrier_groups[[1]]]]) / denominator
        )
    )

save_difference_heatmap(
    diff_df,
    file.path(
        plot_dir,
        sprintf('source_target_heatmap_sig%s_APOE_difference_complexheatmap.pdf', n_samples_sig)
    )
)

session_info()
