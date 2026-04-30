## Louise Huuki-Myers, April 2026
## Build SummarizedExperiment Objetct for Xenium data

#### Set up ####

library("here")
library("SpatialExperiment")
library("SingleCellExperiment")
library("tidyverse")
library("qs2")

data_dir <- here("processed-data", "21_Xenium", "07_xenium_build_spe")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# "/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/raw-data/xenium"
xenium_files <- list.files(here("raw-data", "xenium"))
message("Samples: ", length(xenium_files))

xenium_experiment_table <- tibble(path = xenium_files) |>
    separate(path, into = c(NA, "chip", "sample", "date", "chip2"), sep = "__") |>
    separate(sample, into = c("BrNum", "experiment", "Run"), sep = "_") |>
    mutate(chip = paste0("c", chip)) ## leading 0 protection

xenium_experiment_table |> count(Run, chip)
#   Run   chip         n
# <chr> <chr>    <int>
# 1 Run1  c0102273     4
# 2 Run1  c0102447     4
# 3 Run2  c0061389     4
# 4 Run2  c0102694     4

write_csv(xenium_experiment_table, here(data_dir, "xenium_experiment_details.csv"))

spe_list <- map(xenium_files, function(sample_dir){
    
    message(Sys.time(), sample_dir)
    sample_path = here("raw-data", "xenium", sample_dir)
    
    counts_path <- here(sample_path, "cell_feature_matrix.h5")
    cell_info_path <- here(sample_path, "cells.csv.gz")
    
    spe <- DropletUtils::read10xCounts(counts_path)
    counts(spe) <- methods::as(DelayedArray::realize(counts(spe)), "dgCMatrix") # Convert to delayed array
    
    cell_info <- vroom::vroom(cell_info_path)
    
    colData(spe) <- cbind(colData(spe), cell_info)
    spe <- toSpatialExperiment(spe, spatialCoordsNames = c("x_centroid", "y_centroid"))
    rownames(spe) <- rowData(spe)$Symbol # change rownames to gene symbol
    
    spe$BrNum <- gsub(".*(Br[0-9]{4}).*", "\\1", spe$Sample)
    
    return(spe)
})

#### Combine SPE ####
message(Sys.time(), " - Combine data")
spe <- do.call(cbind, spe_list)

#### dononr info ####
donor_info <- read.csv(here("processed-data", "00_project_prep", "05_pathology", "sample_taupathy.csv")) |>
    tibble::column_to_rownames("BrNum") |>
    mutate(Ancestry = factor(Ancestry),
           APOE = factor(APOE),
           APOE_carrier = factor(APOE_carrier),
           Sex = factor(Sex)
    )


colData(spe) <- cbind(colData(spe), donor_info[spe$BrNum,])

xenium_experiment_table <- xenium_experiment_table |> column_to_rownames("BrNum")

colData(spe) <- cbind(colData(spe), xenium_experiment_table[spe$BrNum,])

## Save data
message(Sys.time(), " - Save qs2")
qs2::qs_save(spe, file = here(data_dir, "spe_xenium.qs2"))


# slurmjobs::job_single('07_xenium_build_spe', create_shell = TRUE, memory = '100G', command = "Rscript 07_xenium_build_spe")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

