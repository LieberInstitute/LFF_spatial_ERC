#   A version of 06_precast.* but using capture area as sample ID and initial
#   array coordinates from spaceranger

library(getopt)
library(sessioninfo)
library(here)
library(PRECAST)
library(HDF5Array)
library(Seurat)
library(tidyverse)
library(Matrix)
library(SpatialExperiment)

# Import command-line parameters
spec <- matrix(
    c(
        "k", "k", "1", "numeric", "Number of clusters",
        "input_genes", "i", "1", "character", "'HVG' or 'SVG'"
    ),
    ncol = 5, byrow = TRUE
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

spe_dir <- here("processed-data", "03_stitching", "spe")
if (opt$input_genes == "HVG") {
    out_path <- here(
        "processed-data", "03_stitching", "precast_out_unstitched",
        sprintf("PRECAST_k%s.csv", opt$k)
    )
    hvg_path = here("processed-data", "03_stitching", "HVGs.txt")
} else {
    out_path <- here(
        "processed-data", "03_stitching", "precast_out_unstitched",
        sprintf("PRECAST_k%s_SVG.csv", opt$k)
    )
    svg_path = here(
        "processed-data", "03_stitching", "nnSVG_out_unstitched",
        "top500_SVGs.txt"
    )
}

set.seed(1)
dir.create(dirname(out_path), showWarnings = FALSE)

spe <- loadHDF5SummarizedExperiment(spe_dir)

#   PRECAST expects array coordinates in 'row' and 'col' columns. Use unstitched
#   coordinates intentionally
spe$row <- spe$array_row_original
spe$col <- spe$array_col_original

#   Create a list of three Seurat objects at the capture-area level
seu_list = lapply(
    unique(spe$capture_area),
    function(capture_area) {
        small_spe = spe[, spe$capture_area == capture_area]
        
        CreateSeuratObject(
            #   Bring into memory to greatly improve speed
            counts = as(assays(small_spe)$counts, "dgCMatrix"),
            meta.data = as.data.frame(colData(small_spe)),
            project = "visiumStitched_brain"
        )
    }
)

#   Run PRECAST using either HVGs or SVGs as input
if (opt$input_genes == "HVG") {
    pre_obj <- CreatePRECASTObject(
        seuList = seu_list,
        selectGenesMethod = NULL,
        customGenelist = readLines(hvg_path),
    )
} else {
    pre_obj <- CreatePRECASTObject(
        seuList = seu_list,
        selectGenesMethod = NULL,
        customGenelist = readLines(svg_path),
        premin.spots = 0,
        postmin.spots = 0
    )
}

pre_obj <- AddAdjList(pre_obj, platform = "Visium")

#   Following https://feiyoung.github.io/PRECAST/articles/PRECAST.BreastCancer.html,
#   which involves overriding some default values, though the implications are not
#   documented
pre_obj <- AddParSetting(
    pre_obj,
    Sigma_equal = FALSE, verbose = TRUE, maxIter = 30
)

#   Fit model
pre_obj <- PRECAST(pre_obj, K = opt$k)
pre_obj <- SelectModel(pre_obj)
pre_obj <- IntegrateSpaData(pre_obj, species = "Human")

#   Extract PRECAST results, clean up column names, and export to CSV
pre_obj@meta.data |>
    rownames_to_column("key") |>
    as_tibble() |>
    select(-orig.ident) |>
    rename_with(~ sub("_PRE_CAST", "", .x)) |>
    write_csv(out_path)

session_info()
