# Sep, 2024 - Louise Huuki-Myers
# Plot uHarmony corrected reduced dim plots for Visium data

library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scater")
library("HDF5Array")
library("viridis")
library("purrr")

plot_dir <- here("plots", "05_spe_correct_cluster", "04_plot_harmony_reduced_dims")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_GLM_Harmony"))
spe

reducedDimNames(spe)
# [1] "10x_pca"          "10x_tsne"         "10x_umap"         "PCA_p1"           "PCA_p2"           "TSNE"            
# [7] "UMAP"             "GLMPCA_approx_1k" "GLMPCA_approx_2k" "GLMPCA_approx_5k" "HARMONY"          "UMAP.HARMONY"    
# [13] "TSNE.HARMONY"   

#### Establish color schemes ####
qc_colors <- c(fold = "#4BA402", out_edge = "#097FE0", vessel ="#BB1EF4", None = "#CCCCCC40")
scran_colors <- c(`TRUE` = "red", `FALSE` = "#CCCCCC40")


## categorical
walk(c("sample_id", "round", "APOE"), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = "UMAP.HARMONY", my_var = .x, sufix = "Harmony"))

walk2(c("qc_anno", "scran_discard"), 
      list(qc_colors, scran_colors), 
      ~my_plot_reduced_dim(spe,prefix = "spe", var_type = "cat", dimred = "UMAP.HARMONY", my_var = .x, sufix = "Harmony", color_pal = .y))
