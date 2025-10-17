library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

num_factors = 7
specific_factor = 'Factor3'
fdr_cutoff = 0.05
z_cutoff_genes = 1.96
max_genes = 2
my_views = c("Oligo.3", "Oligo.4", "Oligo.5", "Excit.L5.2", "WM~Sp09D06")

cell_colors_path = here(
    "processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"
)
spd_colors_path = here("processed-data", "SpD_colors.Rdata")
model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
sce_path = here(
    "processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn",
    "sce_pseudo_DGE-cell_type_anno.RDS"
)
plot_dir = here("plots", "14_MOFA", "04_reduced_gw_heatmap")

dir.create(file.path(plot_dir, 'GO'), recursive = TRUE, showWarnings = FALSE)

#   Prepare view colors
load(cell_colors_path, verbose = TRUE)
load(spd_colors_path, verbose = TRUE)
view_colors = c(cell_type_colors$anno, SpD_colors)[my_views]
view_colors = c(view_colors, Multi = "grey30")

model = load_model(model_path)
sce = readRDS(sce_path)
