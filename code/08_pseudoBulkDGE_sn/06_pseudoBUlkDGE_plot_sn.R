## Louise Huuki-Myers, May 2025
## Plot results of sn pseudoBulkDGE 

library("SingleCellExperiment")
library("scran")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

#### Set up dirs ####
plot_dir <- here("plots", "08_pseudoBulkDGE_sn", "03_pseudoBulkDGE_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS"))


map(de_results, ~min(.x$FDR, na.rm = TRUE))

map_int(de_results, ~sum(!is.na(.x$FDR)))

de_results_tb <- map2_dfr(de.results, names(de.results), ~as.data.frame(.x) |> mutate(cell_type_anno = .y)) |>
    mutate(DE_class = case_when(logFC > 1 & PValue < 0.05 ~ "E4+",
                                logFC < -1 & PValue < 0.05 ~ "E2+",
                                TRUE ~"Other")) |>
    filter(!is.na(logFC))

de_results_tb |> filter(FDR < 0.1)


de_results_tb |> filter(FDR < 0.05) |> dplyr::count(cell_type_anno)
de_results_tb |> filter(PValue < 0.05) |> dplyr::count(cell_type_anno)

de_results_tb |> dplyr::count(cell_type_anno, DE_class)


de_results_tb |> group_by(cell_type_anno) |> dplyr::arrange(FDR) |> dplyr::slice(1)

de_results_tb |> filter(gene_name == "APOE") |> arrange(FDR)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

violin_plot_cell_type_anno <- de_results_tb |>
    ggplot(aes(x = logFC, y = -log10(PValue), color = DE_class)) +
    geom_point(size = .6) +
    geom_text_repel(aes(label = ifelse(PValue < 0.05 | abs(logFC) > 1, gene_name, "")), size = 1.7) +
    facet_wrap(~cell_type_anno) +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed", color = "red") +
    scale_color_manual(values = c(APOE_carrier_colors, Other = "grey")) +
    theme_bw() +
    labs(title = "Visium PseudoBulkDGE", 
         subtitle = as.character(models$carrier)[[2]])

ggsave(violin_plot_cell_type_anno, filename = here(plot_dir, "sn_DGE_violin_plot.png"), width = 12, height = 10)
