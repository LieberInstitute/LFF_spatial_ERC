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

# library("tidyverse")
# library("crumblr")
# library("variancePartition")

data_dir <- here("processed-data", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

## SET UP MULTICORE PARAM
ncores <- 4  # match Sys.getenv('SLURM_CPUS_ON_NODE')
bp <- MulticoreParam(workers = ncores) #or bp <- MulticoreParam(4)


#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

names(rctd_data@results)

head(rctd_data@results$results_df)

## missing cells from ectd data
ncol(spe) - nrow(rctd_data@results$results_df) #17209

setequal(colnames(spe), rownames(rctd_data@results$results_df))

setdiff(rownames(rctd_data@results$results_df), colnames(spe))

spe <- spe[,rownames(rctd_data@results$results_df)]

identical(colnames(spe), rownames(rctd_data@results$results_df))

colData(spe) <- cbind(colData(spe), rctd_data@results$results_df)


#Save SPE with addititional data
message(Sys.time(), " - Saving cleaned SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_cell_types.qs2"))


# get variance explained
pca_var <- attr(reducedDim(spe, "PCA"), "varExplained")
pca_var_pct <- pca_var / sum(pca_var) * 100


#Save cleaned SPE
message(Sys.time(), " - Saving cleaned SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_PCA.qs2"))


# elbow plot
var_explained <- data.frame(PC = seq_along(pca_var_pct),
           var_explained = pca_var_pct) 

var_explained_elbow <- var_explained |>
    ggplot(aes(x = PC, y = var_explained)) +
    geom_point() +
    geom_line() +
    labs(x = "PC", y = "Variance Explained (%)") +
    theme_bw()

ggsave(var_explained_elbow, filename = here(plot_dir, "var_explained_elbow.png"))

source(here("code", "utils", "my_plot_reduced_dim.R"))

#### plot ####
## categorical
walk(c("sample_id", "seq_round", "chip","APOE", "first_type"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat"))
walk(c("sample_id", "seq_round", "chip","APOE", "first_type"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat"))

## continuous
walk(c("sum", "detected", "subsets_Mito_percent"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "con"))

walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", var_type = "express", my_var = .x))
walk(c("MBP", "SNAP25", "SLC17A7" ,"GFAP", "GAD1", "CLDN5", "OLIG2", "TMEM119"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", var_type = "express", my_var = .x))

# slurmjobs::job_single('10_xenium_cell_types', create_shell = TRUE, memory = '100G', command = "Rscript 10_xenium_cell_types.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

