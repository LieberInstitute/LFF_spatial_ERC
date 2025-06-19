## Louise Huuki-Myers, April 2025
## create various plots for the ERC sn data

library("SingleCellExperiment")
# library("jaffelab")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
library("spatialLIBD")
# library("ComplexHeatmap")
# library("bluster")
# library("dendextend")
# library("readxl")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "33_sn_subcluster_summary")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "33_sn_subcluster_summary")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

pd <- as.data.frame(colData(sce))

table(sce$cell_type_anno)

## load colors
# load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors

#### enrichment data ####
erc_sn_modeling_results <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_anno.rds"))

enrichment_top5 <- erc_sn_modeling_results$enrichment

#### Cell Type Summary ####

cell_type_summary <- pd |>
    group_by(cell_type_class, cell_type_broad, cell_type_anno) |>
    summarise(n_nuclei = n(),
              n_donors = length(unique(sample_id)))



#### Quality Metrics ####

sum_umi_violin <- ggplot(pd, aes(x = cell_type_anno, y = sum, fill = cell_type_anno)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = cell_type_colors$anno) +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(sum_umi_violin, filename = here(plot_dir, "ERC_sn_sum_umi_violin.png"), width = 10)

#### violin plots of select genes ####
rownames(sce) <- rowData(sce)$gene_name

apoe_plot <- plot_gene_express(sce, genes = "APOE", 
                               category = "cell_type_anno", 
                               color_pal = cell_type_colors$anno) +
    labs(x = "Cell Type") +
    coord_flip()

ggsave(apoe_plot, filename = here(plot_dir, "ERC_sn_gene_expres_APOE.png"), width = 4)



