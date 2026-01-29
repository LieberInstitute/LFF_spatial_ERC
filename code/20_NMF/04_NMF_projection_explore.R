## Louise Huuki-Myers, Jan 2026
## Correlate NMF factors with variables

#### Set up ####

# library("SpatialExperiment")
library("tidyverse")
library("sessioninfo")
library("here")
library("ComplexHeatmap")
library("spatialLIBD")
library("projectR")

data_dir <- here("processed-data", "20_NMF", "04_NMF_projection_explore")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "20_NMF", "04_NMF_projection_explore")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####

## snRNA-seq
# message(Sys.time(), " - Load HDF5 sce")
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## spatail data
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

## Drop Br1289
# sce <- sce[, sce$BrNum != "Br1289"]
spe <- spe[, spe$BrNum != "Br1289"]

## NMF data 
proj <- readRDS(here("processed-data", "20_NMF", "03_NMF_projection", "ECR_nmf_projection.Rds"))
dim(proj)
proj[1:5,1:5]

ncol(spe) == nrow(proj)

## add Proj to colData
colData(spe) <- cbind(colData(spe),proj)


#### Summarize ####

proj_long <- colData(spe) |>
    as.data.frame() |>
    select(SpD, starts_with("nmf")) |>
    pivot_longer(!SpD, names_to = "NMF", values_to = "proj")

proj_summary <- proj_long |> 
    group_by(SpD, NMF) |> 
    summarise(median = median(proj))

proj_summary |> arrange(-median)

#### Visulize ####

nmf_vis <- vis_gene(spe,
                    geneid = "nmf43")

ggsave(nmf_vis, filename = here(plot_dir, "nmf_vis.png"))


# slurmjobs::job_single('04_NMF_projection_explore', create_shell = TRUE, memory = '200G', command = "Rscript 04_multistudy_Oligo_subcluster.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
