library(here)
library(qs2)
library(sessioninfo)
library(SpatialExperiment)
library(Giotto)

in_path = here(
    'processed-data', '21_Xenium', '10_xenium_cell_types',
    'spe_xenium_cell_types.qs2'
)
out_dir = here('processed-data', '24_xenium_liana', '01_to_python')
python_path = '/users/neagles/.conda/envs/giotto_env/bin/python'

GiottoClass::set_giotto_python_path(python_path)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spe = qs_read(in_path)

#   Trim down the object
assays(spe) = list(logcounts = assays(spe)$logcounts)

spe |>
    spatialExperimentToGiotto() |>
    giottoToAnnData(save_directory = out_dir)

session_info()
