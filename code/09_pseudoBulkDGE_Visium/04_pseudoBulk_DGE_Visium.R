## Louise Huuki-Myers, April 2025
## Run  pseudoBulkDGE

library("spatialExperiment")
library("edgeR")
library("scran")
library("tidyverse")
library("here")
library("sessioninfo")

#### Set up dirs ####
plot_dir <- here("plots", "09_pseudoBulkDGE_Visium", "03_pseudoBulk_DGE_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
spe_pb <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))


model.matrix(~APOE_carrier + Sex + Age, colData(spe_pb))

de_results <- pseudoBulkDGE(
    spe_pb,
    label = spe_pb$SpD,
    design = ~APOE_carrier + Anc_Afr + Age + Sex + Rin + Visium_slide,
    coef = "APOE_carrierE4+",
    condition = "APOE_carrier"
    row.data = rowData(spe_pb),
    method = "edgeR"
)

map(de_results, ~sum(.x$FDR<0.05, na.rm = TRUE))
map(de_results, ~sum(.x$FDR<0.2, na.rm = TRUE))

map(de_results, ~min(.x$FDR, na.rm = TRUE))

map_int(de_results, ~sum(!is.na(.x$FDR)))

de_results_tb <- map2_dfr(de_results, names(de_results), ~as.data.frame(.x) |> mutate(SpD = .y)) |>
    mutate(DE_class = case_when(logFC > 1 & PValue < 0.05 ~ "E4+",
                                logFC < -1 & PValue < 0.05 ~ "E2+",
                                TRUE ~"Other")) |>
    filter(!is.na(logFC))

de_results_tb |> filter(FDR < 0.05) |> dplyr::count(SpD)
de_results_tb |> filter(PValue < 0.05) |> dplyr::count(SpD)

de_results_tb |> dplyr::count(SpD, DE_class)


de_results_tb |> group_by(SpD) |> dplyr::arrange(FDR) |> dplyr::slice(1)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

(pval_lim <- min(de_results_tb$PValue[de_results_tb$FDR < 0.05]))
    
violin_plot_SpD <- de_results_tb |>
    ggplot(aes(x = logFC, y = -log10(PValue), color = DE_class)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = ifelse(PValue < 0.05 | abs(logFC) > 1, gene_name, "")), size = 2) +
    facet_wrap(~SpD) +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed", color = "red") +
    scale_color_manual(values = c(APOE_carrier_colors, Other = "grey")) +
    theme_bw() +
    labs(title = "Visium PseudoBulkDGE", 
         subtitle = "~APOE_carrier + Anc_Afr + Age + Sex + Rin + Visium_slide")

ggsave(violin_plot_SpD, filename = here(plot_dir, "Visium_DGE_violin_plot.png"), width = 10, height = 8)
    
