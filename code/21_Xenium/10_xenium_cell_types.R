## Louise Huuki-Myers, May 2026
## Explore RCTD results

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("BiocParallel")
library("scran")
library("scater")

library("tidyverse")
# library("crumblr")
# library("variancePartition")

data_dir <- here("processed-data", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

## SET UP MULTICORE PARAM
# ncores <- 4  # match Sys.getenv('SLURM_CPUS_ON_NODE')
# bp <- MulticoreParam(workers = ncores) #or bp <- MulticoreParam(4)
# 

#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))


#### Add RCTD cell types to spe ####
message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

names(rctd_data@results)

head(rctd_data@results$results_df)

## missing cells from rctd data
ncol(spe) - nrow(rctd_data@results$results_df) #17209

spe$cell_id <- colnames(spe)

rctd_tb <- rctd_data@results$results_df |>
    tibble::rownames_to_column("cell_id")

# extract colData with cell IDs
col_df <- as.data.frame(colData(spe)) |>
    left_join(rctd_tb, by = "cell_id")
    
# put back into spe
rownames(col_df) <- col_df$cell_id
colData(spe) <- DataFrame(col_df)


table(is.na(spe$spot_class), spe$sum_gex < 100)

#         FALSE   TRUE
# FALSE 450307      0
# TRUE       0  17209

## 3% of nuclei <100 counts (sum_gex) and missing RCTD nuclei 
sum(is.na(spe$spot_class))/ncol(spe)

#Save SPE with additional data
message(Sys.time(), " - Saving cleaned SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_cell_types.qs2"))

# spe <- qs_read(here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_PCA.qs2"))

source(here("code", "utils", "my_plot_reduced_dim.R"))

#### plot ####
## categorical
walk2(c("first_type"), list(cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat", color_pal = .y))
walk2(c("first_type"), list(cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", color_pal = .y))

walk2(c("second_type"), list(cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", color_pal = .y))

walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat"))
walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat"))

walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat"))
walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat"))





# slurmjobs::job_single('10_xenium_cell_types', create_shell = TRUE, memory = '100G', command = "Rscript 10_xenium_cell_types.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

