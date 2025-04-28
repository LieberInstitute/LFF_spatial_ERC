## Louise Huuki-Myers, April 2025
## create various plots for the ERC sn data

library("SingleCellExperiment")
# library("jaffelab")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
# library("ComplexHeatmap")
# library("bluster")
# library("dendextend")
# library("readxl")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "28_sn_clean_plot")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "28_sn_clean_plot")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

pd <- as.data.frame(colData(sce))

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Quality Metrics ####

sum_umi_violin <- ggplot(pd, aes(x = cell_type_anno, y = sum, fill = cell_type_anno)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = cell_type_colors$anno) +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(sum_umi_violin, filename = here(plot_dir, "ERC_sn_sum_umi_violin.png"), width = 10)


