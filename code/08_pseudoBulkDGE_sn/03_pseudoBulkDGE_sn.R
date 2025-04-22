## Louise Huuki-Myers, April 2025
## Run  pseudoBulkDGE on cell types from sn data

library("SingleCellExperiment")
# library("edgeR")
library("scran")
library("tidyverse")
library("here")
library("sessioninfo")

#### Set up dirs ####
plot_dir <- here("plots", "08_pseudoBulkDGE_sn", "03_pseudoBulkDGE_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data","sce_pseudo_DGE.RDS"))
dim(sce_pb)

sce_pb <- sce_pb[,!sce_pb$cell_type_anno %in% c("Inhib.Pax6", "Excit.L5_6_NP")]
sce_pb$cell_type_anno <- droplevels(sce_pb$cell_type_anno)

table(sce_pb$cell_type_anno)

#### Design model + run DGE ####
model.matrix(~APOE_carrier + Anc_Afr + Age + Sex + Rin + exp_round, colData(sce_pb))

de_results <- pseudoBulkDGE(
    sce_pb,
    label = sce_pb$cell_type_anno,
    design = ~APOE_carrier + Anc_Afr + Age + Sex + exp_round,
    coef = "APOE_carrierE4+",
    row.data = rowData(sce_pb),
    method = "edgeR"
)

# error when including Rin ?
# Error in glmFit.default(sely, design, offset = seloffset, dispersion = 0.05,  : 
#                             nrow(design) disagrees with ncol(y)

map(de_results, ~sum(.x$FDR<0.05, na.rm = TRUE))
map(de_results, ~sum(.x$FDR<0.2, na.rm = TRUE))

map(de_results, ~min(.x$FDR, na.rm = TRUE))

map_int(de_results, ~sum(!is.na(.x$FDR)))

de_results_tb <- map2_dfr(de_results, names(de_results), ~as.data.frame(.x) |> mutate(cell_type_anno = .y)) |>
    mutate(DE_class = case_when(logFC > 1 & PValue < 0.05 ~ "E4+",
                                logFC < -1 & PValue < 0.05 ~ "E2+",
                                TRUE ~"Other")) |>
    filter(!is.na(logFC))

de_results_tb |> filter(FDR < 0.05) |> dplyr::count(cell_type_anno)
de_results_tb |> filter(PValue < 0.05) |> dplyr::count(cell_type_anno)

de_results_tb |> dplyr::count(cell_type_anno, DE_class)


de_results_tb |> group_by(cell_type_anno) |> dplyr::arrange(FDR) |> dplyr::slice(1)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

(pval_lim <- min(de_results_tb$PValue[de_results_tb$FDR < 0.05]))

violin_plot_cell_type_anno <- de_results_tb |>
    ggplot(aes(x = logFC, y = -log10(PValue), color = DE_class)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = ifelse(PValue < 0.05 | abs(logFC) > 1, gene_name, "")), size = 2) +
    facet_wrap(~cell_type_anno) +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed", color = "red") +
    scale_color_manual(values = c(APOE_carrier_colors, Other = "grey")) +
    theme_bw() +
    labs(title = "Visium PseudoBulkDGE", 
         subtitle = "~APOE_carrier + Anc_Afr + Age + Sex + Rin + Visium_slide")

ggsave(violin_plot_cell_type_anno, filename = here(plot_dir, "Visium_DGE_violin_plot.png"), width = 10, height = 8)

