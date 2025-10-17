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
highlighted_views = factor_weights |>
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
    slice_head(n = num_views) |>
    pull(view)

stopifnot(all(highlighted_views %in% colnames(factor_weights)))
factor_weights = factor_weights[, highlighted_views]

col_fun_r2 = colorRamp2(
    seq(0, min(100, max(factor_weights, na.rm = TRUE) + 5), length = 50),
    hcl.colors(50, "Oranges", rev = TRUE)
)

################################################################################
#   Donor-level covariates column annotation
################################################################################

#   Test association of a couple donor-level covariates with each factor
assoc_list = list()
for (var_name in c("APOE_carrier", "taupathy")) {
    assoc_list[[var_name]] = get_associations(
        model = model,
        metadata = samples_metadata(model),
        sample_id_column = "sample",
        test_variable = var_name,
        test_type = "categorical",
        group = FALSE
    )
}

#   Show unadjusted p-value for consistency with other plots
for (var_name in names(assoc_list)) {
    assoc_list[[var_name]]$adj_pvalue = assoc_list[[var_name]]$p.value
}

assoc_pvals = assoc_list |>
    tibble::enframe(name = "test") |>
    tidyr::unnest(c(value)) |>
    dplyr::mutate(log_adjpval = -log10(.data$adj_pvalue)) |>
    dplyr::select(test, Factor, log_adjpval) |>
    tidyr::pivot_wider(names_from = test, values_from = log_adjpval) |>
    dplyr::select(-Factor) |>
    as.matrix()

col_fun_assoc = colorRamp2(
    seq(0, max(assoc_pvals) + 0.5, length = 20),
    hcl.colors(20, "Purples", rev = TRUE)
)

column_ha = HeatmapAnnotation(
    "R2" = factor_weights,
    "pvalue" = assoc_pvals,
    gap = unit(2.5, "mm"),
    border = TRUE,
    col = list(R2 = col_fun_r2, pvalue = col_fun_assoc)
)

################################################################################
#   Main heatmap body
################################################################################

factor_matrix = get_factors(model, factors = "all")$single_group

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

pdf(plot_path, width = 4, height = as.integer(round(num_views * 0.75) + 1))
ComplexHeatmap::draw(scores_hmap, padding = unit(c(2, 2, 2, 15), "mm"))
dev.off()

session_info()
