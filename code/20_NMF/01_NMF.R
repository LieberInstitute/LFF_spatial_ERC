## Louise Huuki-Myers, Jan 2026
## Run NMF on ERC snRNA-seq data
## following https://github.com/LieberInstitute/spatialdACC/blob/main/code/snRNA-seq/06_NMF/01_NMF_factor.R

#### Set up ####

library("singlet")
library("Seurat")
library("dplyr")
library("ggplot2")
library("sessioninfo")
library("here")

data_dir <- here("processed-data", "01_NMF", "01_NMF")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "01_NMF", "01_NMF")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)



#### load data ####
pbmc3k <- singlet::get_pbmc3k_data()

pbmc3k
class(pbmc3k)

set.seed(123)
pbmc3k <- NormalizeData(pbmc3k) |>
    singlet::RunNMF()
