## October 2024, Louise Huuki-Myers
## Run PRECAST clustering
## Adapted from https://github.com/LieberInstitute/visiumStitched_brain/blob/9d2f716bd76a3359b2aeca1cac3c90720c12bfb1/code/03_stitching/07_precast_unstitched.R

library("getopt")
library("sessioninfo")
library("here")
library("PRECAST")
library("HDF5Array")
library("Seurat")
library("tidyverse")
library("Matrix")
library("SpatialExperiment")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "06_PRECAST")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
# 
# plot_dir <- here("plots", "05_spe_correct_cluster", "06_PRECAST")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
spec <- matrix(
    c(
        "k", "k", "1", "numeric", "Number of clusters",
        "input_genes", "i", "1", "character", "'HVG' or 'SVG'"
    ),
    ncol = 5, byrow = TRUE
)
opt <- getopt(spec)

#test
opt <- list(k=2, input_genes = "HVG")

print("Using the following parameters:")
print(opt)

## load HVGs (TODO SVGs)
load(here("processed-data", "02_build_spe","01_preprocess_spe", "top_hvgs.Rdata"), verbose = TRUE)
# top.hvgs
# map_int(top.hvgs, length)
# p1   p2 
# 2021 4042

## define output path
out_path <- here(
    data_dir,
    sprintf("PRECAST_%s_k%02d.csv", opt$input_genes, opt$k)
)

dir.create(dirname(out_path), showWarnings = FALSE)

#### Load SPE data ####
message(Sys.time(), " - Loading SPE data")
spe <- loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))

#   PRECAST expects array coordinates in 'row' and 'col' columns. Use unstitched
#   coordinates intentionally
spe$row <- spe$array_row
spe$col <- spe$array_col

#### convert to seurat object ####
#   Create a list of three Seurat objects at the capture-area level
message(Sys.time(), " - Convert to Seurat Object")
seu_list = lapply(
    unique(spe$sample_id),
    function(sample_id) {
        small_spe = spe[, spe$sample_id == sample_id]
        
        CreateSeuratObject(
            #   Bring into memory to greatly improve speed
            counts = as(assays(small_spe)$counts, "dgCMatrix"),
            meta.data = as.data.frame(colData(small_spe)),
            project = "LFF_ERC"
        )
    }
)

####   Run PRECAST using either HVGs or SVGs as input ####
set.seed(1)
message(Sys.time(), " - Create PRECAST object")
if (opt$input_genes == "HVG") {
    pre_obj <- CreatePRECASTObject(
        seuList = seu_list,
        selectGenesMethod = NULL,
        customGenelist = top.hvgs$p1,
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
message(Sys.time(), " - Run PRECAST, k=", opt$k)
pre_obj <- PRECAST(pre_obj, K = opt$k)
pre_obj <- SelectModel(pre_obj)
pre_obj <- IntegrateSpaData(pre_obj, species = "Human")

#   Extract PRECAST results, clean up column names, and export to CSV
message(Sys.time(), " - Extract & export results")
pre_obj@meta.data |>
    rownames_to_column("key") |>
    as_tibble() |>
    select(-orig.ident) |>
    rename_with(~ sub("_PRE_CAST", "", .x)) |>
    write_csv(out_path)

# slurmjobs::job_loop(loops = list(input_genes = c("HVG")),
#                     name ='06_PRECAST', create_shell = TRUE, memory = '50G')

# slurmjobs::array_submit(name = "06_PRECAST", task_ids = c(3:10), submit = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
