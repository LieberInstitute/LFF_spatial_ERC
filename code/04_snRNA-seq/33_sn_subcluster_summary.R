## Louise Huuki-Myers, April 2025
## create various plots for the ERC sn data

library("SingleCellExperiment")
library("jaffelab")
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

dim(sce)

pd <- as.data.frame(colData(sce))

table(sce$cell_type_anno)

## load colors
# load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
cell_type_colors <- metadata(sce)$cell_type_colors

#### enrichment data ####
## global
enrichment_stats_top <- read.csv(here("processed-data", "04_snRNA-seq", "31_sn_subcluster_heatmap", "sce_subcluster_enrichment_top5.csv"))

enrichment_stats_top_list <- enrichment_stats_top |>
    group_by(cell_type_anno) |>
    summarise(top_enrichment_genes_global = paste(ensembl, collapse = ", "))

#### Cell Type Summary ####
cell_type_summary <- pd |>
    group_by(cell_type_class, cell_type_broad, cell_type_anno) |>
    summarise(n_nuclei = n(),
              prop = n_nuclei/ncol(sce),
              n_donors = length(unique(sample_id)),
              median_sum = median(sum),
              median_detected = median(detected),
              median_Mito_percent = median(subsets_Mito_percent),
              median_scDblFinder.score = median(scDblFinder.score)) |>
    left_join(enrichment_stats_top_list)

cell_type_summary |> filter(cell_type_anno == "Inhib.Pax6")
# cell_type_summary |> filter(cell_type_anno == "Oligo.3")
# cell_type_summary |> filter(cell_type_anno == "Astro.2")

write.csv(cell_type_summary, file = here(data_dir, "ERC_sn_subcluster_summary_cell_type.csv"))

#### Donor Summary ####

sample_summary <- pd |>
    group_by(sample_id, exp_round, seq_round) |>
    summarise(n_nuclei = n(),
              median_sum = median(sum),
              median_detected = median(detected),
              median_Mito_percent = median(subsets_Mito_percent),
              median_scDblFinder.score = median(scDblFinder.score)) 

summary(sample_summary)

write.csv(sample_summary, file = here(data_dir, "ERC_sn_subcluster_summary_sample.csv"))

#### Quality Metrics ####

qc_violin_sum_umi <- ggplot(pd, aes(x = cell_type_anno, y = sum, fill = cell_type_anno)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = cell_type_colors$anno) +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    facet_grid(.~cell_type_broad, scales = "free_x", space = "free") +
    labs(x = "Cell Type", y = "log10(Sum UMI)") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(qc_violin_sum_umi, filename = here(plot_dir, "ERC_sn_subcluster_QC_violin_sum_umi.png"), width = 7, height =4)

qc_violin_detected <- ggplot(pd, aes(x = cell_type_anno, y = detected, fill = cell_type_anno)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = cell_type_colors$anno) +
    # scale_y_continuous(trans='log10') +
    theme_bw() +
    facet_grid(.~cell_type_broad, scales = "free_x", space = "free") +
    labs(x = "Cell Type", y = "Detected Genes") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(qc_violin_detected, filename = here(plot_dir, "ERC_sn_subcluster_QC_violin_detected.png"), width = 7, height =4)


qc_violin_mito <- ggplot(pd, aes(x = cell_type_anno, y = subsets_Mito_percent, fill = cell_type_anno)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = cell_type_colors$anno) +
    # scale_y_continuous(trans='log10') +
    theme_bw() +
    facet_grid(.~cell_type_broad, scales = "free_x", space = "free") +
    labs(x = "Cell Type", y = "Precent Mito") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(qc_violin_mito, filename = here(plot_dir, "ERC_sn_subcluster_QC_violin_Mito_percent.png"), width = 7, height =4)

#### n cell barplot ####
n_cell_barplot <- cell_type_summary |>
    ggplot(aes(x = cell_type_anno, y = n_nuclei, fill = cell_type_anno)) +
    geom_col() +
    geom_text(aes(label = n_nuclei), size = 2, vjust = -0.5) +
    scale_fill_manual(values = cell_type_colors$anno, guide = "none") +
    facet_grid(.~cell_type_broad, scales = "free_x", space="free_x") +
    labs(x = "Cell Type", y = "n Nuclei") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(n_cell_barplot, filename = here(plot_dir, "ERC_sn_subcluster_n_nuclei_barplot.png"), width = 9, height = 5)


#### violin plots of select genes ####
rownames(sce) <- rowData(sce)$gene_name

plot_one_gene <- function(gene){
    
    apoe_plot <- plot_gene_express(sce, genes = gene, 
                                   category = "cell_type_anno", 
                                   color_pal = cell_type_colors$anno) +
        labs(x = "Cell Type") 
    
    ggsave(apoe_plot, filename = here(plot_dir, sprintf("ERC_sn_gene_expres_%s.png", gene), height = 4))
    
}

apoe_plot <- plot_gene_express(sce, genes = "APOE", 
                               category = "cell_type_anno", 
                               color_pal = cell_type_colors$anno) +
    labs(x = "Cell Type") +
    coord_flip()

ggsave(apoe_plot, filename = here(plot_dir, "ERC_sn_gene_expres_APOE.png"), width = 4)

#### plot marker genes ####
message(Sys.time(), " - Plot marker genes")

## read in marker genes from lit
lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    mutate(in_data = gene_name %in% rowData(sce)$gene_name,
           cell_type_broad = factor(cell_type_broad, levels = levels(sce$cell_type_broad))) |>
    arrange(cell_type_broad)

lit_markers <- lit_markers |> filter(in_data)

lit_markers_list <- map(splitit(lit_markers$cell_type), ~lit_markers$gene_name[.x])

lit_marker_order <- lit_markers |> count(cell_type_broad, cell_type) |> pull(cell_type)
lit_markers_list <- lit_markers_list[lit_marker_order]

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_sn_subtype_lit_markers.pdf"),
                         cellType_col = "cell_type_anno",
                         gene_name_col = "gene_name",
                         color_pal = cell_type_colors$anno,
)


