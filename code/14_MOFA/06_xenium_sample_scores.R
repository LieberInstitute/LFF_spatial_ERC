library(SpatialExperiment)
library(MOFAcellulaR)
library(MOFA2)
library(tidyverse)
library(here)
library(sessioninfo)

spe_path = here(
    'processed-data', '21_Xenium', '14_xenium_pseudobulk_model_register',
    'spe_xenium_pseudobulk-cell_type_anno.rds'
)
model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
sn_factor_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine",
    "MOFA_factor_df-sn_fine.csv"
)
out_dir = here(
    "processed-data", "14_MOFA", "06_xenium_sample_scores"
)

dir.create(out_dir, showWarnings = FALSE)

spe = readRDS(spe_path)
model = load_model(model_path)

xenium_dat = create_init_exp(
        counts = assays(spe)$counts, coldata = as.data.frame(colData(spe))
    ) |>
    filt_profiles(
        #   Retain all cell types initially
        cts = unique(spe$cell_type_anno),
        ncells = 0,
        counts_col = "ncells",
        ct_col = "cell_type_anno"
    ) |>
    filt_views_bysamples(nsamples = 2) |> #   Drop cell types seen in very few samples
    filt_gex_byexpr(min.count = 5, min.prop = 0.25) |> #   Require a gene to have 5 counts in at least 25% of samples
    filt_views_bygenes(ngenes = 15) |> #   Drop cell types having below 15 genes at this point
    filt_samples_bycov(prop_coverage = 0.9) |> #   Drop samples below 90% coverage of genes
    tmm_trns(scale_factor = 1000000) |> #   Normalize expression
    filt_views_bygenes(ngenes = 15) |>  #   Again, drop cell types having below 15 genes at this point
    center_views() |> # recommended specifically when projecting to new data
    pb_dat2MOFA(sample_column = 'BrNum')

projected_df = project_data(model = model, test_data = xenium_dat) |>
    as.data.frame() |>
    rownames_to_column('donor') |>
    as_tibble()

write_csv(projected_df, file.path(out_dir, "xenium_factor_scores.csv"))

sn_df = read_csv(sn_factor_path, show_col_types = FALSE) |>
    filter(Factor == 'Factor3') |>
    select(sample, value) |>
    dplyr::rename(donor = sample, sn_factor_score = value)

projected_df |>
    select(donor, Factor3) |>
    dplyr::rename(xen_factor_score = Factor3) |>
    inner_join(sn_df, by = 'donor') |>