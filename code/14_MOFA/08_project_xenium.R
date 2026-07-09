library(here)
library(tidyverse)
library(sessioninfo)
library(qs2)
library(SpatialExperiment)
library(spatialLIBD)
library(MOFAcellulaR)
library(MOFA2)

spe_path = here(
    'processed-data', '21_Xenium', '13_xenium_bansky_embedding',
    'spe_xenium_bansky.qs2'
)
astro_path = here(
    'processed-data', '21_Xenium', '20_xenium_Oligo3_Astro',
    'Oligo3_Astro_neighbor_df.Rdata'
)
task_path = here(
    'processed-data', '14_MOFA', '07_build_tasks', 'task_map.csv'
)
model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
out_dir = here("processed-data", "14_MOFA", "08_project_xenium")

this_task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

dir.create(out_dir, showWarnings = FALSE)

################################################################################
#   Load, filter, and pseudobulk data according to the array task
################################################################################

spe = qs_read(spe_path)

astro_df = load(astro_path) |>
    get() |>
    dplyr::rename(cell_id = reference_barcode) |>
    select(cell_id, APOE_level)

spe$astro_group = tibble(cell_id = spe$cell_id) |>
    left_join(astro_df, by = 'cell_id') |>
    pull(APOE_level)

task_df = read_csv(task_path, show_col_types = FALSE) |>
    filter(task_id == this_task_id)

#   Filter Oligo.3 as applicable
if (task_df$astro_group != 'all') {
    spe = spe[
        ,
        (spe$astro_group == task_df$astro_group) &
        (spe$cell_type_anno != 'Oligo.3')
    ]
}

#   Filter spatial domain as applicable
if (task_df$domain != 'all') {
    spe = spe[, spe$SpX == task_df$domain]
}

spe_pb = registration_pseudobulk(
    spe, var_registration = "cell_type_anno", var_sample_id = "sample_id"
)
rownames(spe_pb) = rowData(spe_pb)$gene_id

################################################################################
#   Project MOFA factors to pseudobulked Xenium data
################################################################################

xenium_dat = create_init_exp(
        counts = assays(spe_pb)$counts, coldata = as.data.frame(colData(spe_pb))
    ) |>
    filt_profiles(
        #   Retain all cell types initially
        cts = unique(spe_pb$cell_type_anno),
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
    pb_dat2MOFA(sample_column = 'sample_id')

model = load_model(model_path)

projected_df = project_data(model = model, test_data = xenium_dat) |>
    as.data.frame() |>
    rownames_to_column('donor') |>
    as_tibble() |>
    dplyr::rename(xen_factor_score = Factor3) |>
    select(donor, xen_factor_score)

################################################################################
#   Merge with sn/Visium scores + metadata, and export data for plotting
################################################################################

factor_df = get_tidy_factors(
        model = model, metadata = samples_metadata(model), factor = 'Factor3',
        sample_id_column = "sample"
    ) |>
    dplyr::rename(donor = sample, sn_factor_score = value) |>
    select(donor, sn_factor_score, APOE_carrier, taupathy) |>
    inner_join(projected_df, by = 'donor') |>
    mutate(domain = task_df$domain, astro_group = task_df$astro_group) |>
    write_csv(
        file.path(out_dir, sprintf('factor_data_task%d.csv', this_task_id))
    )

session_info()
