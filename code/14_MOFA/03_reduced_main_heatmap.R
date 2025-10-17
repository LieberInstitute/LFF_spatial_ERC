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