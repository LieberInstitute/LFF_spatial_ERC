## Louise Huuki-Myers, Dec 2025
## Plot gene expression of other Oligo clusters

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
# library("scDotPlot")

data_dir <- here("processed-data", "19_other_Oligo", "03.5_other_Oligo_compile")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "03.5_other_Oligo_compile")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


cor_layer_fn <- list.files( here("processed-data", "19_other_Oligo", "02_other_Oligo_model"), pattern = "cor_layer_anno")

cor_layer_fn <-  c(
    # "other_Oligo_cor_layer_anno_spatialdACC_wOPC_k20.Rds",
    # "other_Oligo_cor_layer_anno_spatialDLPFC_k10.Rds",
    # "other_Oligo_cor_layer_anno_spatialDLPFC_k30.Rds",
    # "other_Oligo_cor_layer_anno_spatialDLPFC_wOPC_k10_OPC.Rds",
    # "other_Oligo_cor_layer_anno_spatialDLPFC_wOPC_k10.Rds",
    # "other_Oligo_cor_layer_anno_spatialHPC_k10.Rds",
   HPC = "other_Oligo_cor_layer_anno_spatialHPC_k20.Rds",
   dACC = "other_OligoOPC_cor_layer_anno_spatialdACC_wOPC_k20.Rds",
   DLPFC = "other_OligoOPC_cor_layer_anno_spatialDLPFC_wOPC_k10.Rds")

names(cor_layer_fn) <- gsub("other_.*?_anno_|.Rds", "", cor_layer_fn)

cor_layer_data <- map(cor_layer_fn, ~readRDS(here("processed-data", "19_other_Oligo", "02_other_Oligo_model", .x)))

pdf(here(plot_dir, "layer_stat_cor_plot_all.pdf"))
walk(cor_layer_data, ~print(layer_stat_cor_plot(cor_stats_layer = .x$cor_layer, annotation = .x$anno)))
dev.off()


         