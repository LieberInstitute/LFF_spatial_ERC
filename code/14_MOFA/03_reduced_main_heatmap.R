#   Produce a custom heatmap showing a subset of views in the column annotation
#   and just APOE_carrier and taupathy covariates

library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
plot_path = here("plots", "14_MOFA", "03_reduced_main_heatmap", "heatmap.pdf")
out_path = here(
    "processed-data", "14_MOFA", "03_reduced_main_heatmap",
    "variance_explained.csv"
)
num_views = 5
specific_factor = "Factor3"

dir.create(dirname(plot_path), showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   R2 column annotation
################################################################################

model = load_model(model_path)
factor_weights = model@cache$variance_explained$r2_per_factor$single_group

factor_weights |>
    as.data.frame() |>
    rownames_to_column('factor_num') |>
    write_csv(out_path)

#   For the specific factor of interest, show just the top few views by R^2
highlighted_views_df = factor_weights |>
    as.data.frame() |>
    rownames_to_column('factor_num') |>
    as_tibble() |>
    pivot_longer(
        cols = -factor_num,
        names_to = 'view',
        values_to = 'r2'
    ) |>
    filter(factor_num == specific_factor) |>
    arrange(desc(r2)) |>
    slice_head(n = num_views) 

highlighted_views <- highlighted_views_df|>
    pull(view)

stopifnot(all(highlighted_views %in% colnames(factor_weights)))
factor_weights = factor_weights[
    specific_factor, highlighted_views, drop = FALSE
]

col_fun_r2 = colorRamp2(
    seq(0, min(100, max(factor_weights, na.rm = TRUE)), length = 50),
    hcl.colors(50, "Oranges", rev = TRUE)
)

column_ha = HeatmapAnnotation(
    "R2" = factor_weights,
    gap = unit(2.5, "mm"),
    border = TRUE,
    col = list(R2 = col_fun_r2)
)

##
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

names(SpD_colors) = sub('~', '_', names(SpD_colors))
view_colors = c(cell_type_colors$anno, SpD_colors)[highlighted_views]

highlighted_views_col_plot <- highlighted_views_df |>
    mutate(view = factor(view, levels = highlighted_views)) |>
    ggplot(aes(x = view, y = r2, fill = view)) +
    geom_col() +
    scale_fill_manual(values = view_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(y = "Factor 3 r2")

ggsave(highlighted_views_col_plot, filename = here("plots", "14_MOFA", "03_reduced_main_heatmap", "factor3_highlighted_views_col_plot.pdf"), height = 2.5, width = 3)

################################################################################
#   Main heatmap body
################################################################################

factor_matrix = get_factors(model, factors = "all")$single_group[
    , specific_factor, drop = FALSE
]

max_fact = abs(max(factor_matrix))
col_fun_fact = colorRamp2(
    seq((max_fact + 0.5) * -1, max_fact + 0.5, length = 50),
    hcl.colors(50, "Green-Brown", rev = TRUE)
)

################################################################################
#   Putting it all together
################################################################################

scores_hmap = Heatmap(
    factor_matrix,
    name = "factor_scores",
    top_annotation = column_ha,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    border = TRUE,
    gap = unit(2.5, "mm"),
    col = col_fun_fact
)

pdf(plot_path, width = 2, height = 3)
ComplexHeatmap::draw(scores_hmap, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

session_info()
