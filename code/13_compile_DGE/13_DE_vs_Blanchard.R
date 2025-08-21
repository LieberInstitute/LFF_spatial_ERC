## Louise Huuki-Myers, Aug 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("jaffelab")
library("ggrepel")

data_dir <- here("processed-data", "13_compile_DGE", "13_DE_vs_Blanchard")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "13_DE_vs_Blanchard")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## sn fine 
DE_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds"))

DE_Oligo <- DE_data |> filter(grepl("Oligo", cluster))


## Blanchard data 
suppressMessages(source(here("external-data", "Blanchard2022", "get_Blanchard_data.R")))

# Blancard_DE_ipsc_Oligo
# Blancard_DE_pm_Oligo
# Blancard_DE_Nebula



DE_Oligo <- DE_data |> filter(grepl("Oligo", cluster))



compare_stats_scatter(carrier_data, mX= "dream", mY="pbDGE", model_name = "carrier")


#### combine data ####
Blancard_DE_pm_Oligo  <- Blancard_DE_pm_Oligo |> select(Blanchard_test = grp2, gene_name = gene, fdr_Blanchard = padj, logFC_Blanchard = logFC)

dge_combined <- DE_Oligo  |>
    select(cluster, contrast, gene_id, gene_name, vlmf_t, fdr_ERC = vlmf_adj.P.Val, logFC_ERC = vlmf_logFC) |>
    inner_join(Blancard_DE_pm_Oligo, by = join_by(gene_name), relationship = "many-to-many")

## calc correlation
cor <- dge_combined |> 
    filter(!is.na(fdr_ERC) & !is.na(fdr_Blanchard)) |>
    group_by(cluster, Blanchard_test) |>
    summarise(cor = cor.test(fdr_ERC, fdr_Blanchard)$estimate,
              p_val = cor.test(fdr_ERC, fdr_Blanchard)$p.value) |>
    ungroup() |>
    mutate(p_val_adj = p.adjust(p_val, method = "bonf"),
           signif = ifelse(p_val_adj < 0.05, "*", ""),
           anno = sprintf("cor=%.2f%s", cor, signif))

cor |> filter(cluster == "Oligo.3")

    
## create scatter plot
source(here("code", "13_compile_DGE", "compare_stats_scatter.R"))

stat_scatter <- compare_stats_scatter(dge_tb = dge_combined, 
                                      stat = "logFC", 
                                      mX = "ERC", 
                                      mY = "Blanchard", 
                                      FDR_cut_mX = 0.05, 
                                      FDR_cut_mY = 0.05, 
                                      model_name = "Blancard PM") +
    facet_grid(cluster~Blanchard_test) +
    geom_label(
        data = cor, ggplot2::aes(x = -Inf, y = Inf, label = anno),
        color = "black",
        alpha = 0.5,
        vjust = "inward", 
        hjust = "inward", 
        size = 2.5
    )

ggsave(stat_scatter, filename = here(plot_dir, "Blanchard_pm_logFC_scatter_Oligo.png"), height = 10, width = 10)


stat_scatter_Oligo.3 <- compare_stats_scatter(dge_tb = dge_combined |> filter(cluster == "Oligo.3") |> mutate(cluster = droplevels(cluster)), 
                                      stat = "logFC", 
                                      mX = "ERC", 
                                      mY = "Blanchard", 
                                      FDR_cut_mX = 0.05, 
                                      FDR_cut_mY = 0.05, 
                                      model_name = "Blancard PM") +
    facet_grid(cluster~Blanchard_test) +
    geom_label(
        data = cor |> filter(cluster == "Oligo.3"), 
        ggplot2::aes(x = -Inf, y = Inf, label = anno),
        color = "black",
        alpha = 0.5,
        vjust = "inward", 
        hjust = "inward", 
        size = 2.5
    )
    
ggsave(stat_scatter_Oligo.3, filename = here(plot_dir, "Blanchard_pm_logFC_scatter_Oligo.3.png"), height = 6, width = 10)

## filter to genes of interest 
suppressMessages(source(here("external-data", "Blanchard2022", "get_Blanchard_gene_list.R")))
# Blanchard_gene_list

Blanchard_genes <- unlist(Blanchard_gene_list)

stat_scatter_Oligo.3_select <- compare_stats_scatter(dge_tb = dge_combined |> 
                                                  filter(cluster == "Oligo.3",
                                                         gene_name %in% Blanchard_genes) |> 
                                                  mutate(cluster = droplevels(cluster)), 
                                      stat = "logFC", 
                                      mX = "ERC", 
                                      mY = "Blanchard", 
                                      FDR_cut_mX = 0.05, 
                                      FDR_cut_mY = 0.05, 
                                      model_name = "Blancard PM") +
    facet_grid(cluster~Blanchard_test) 


ggsave(stat_scatter_Oligo.3_select, filename = here(plot_dir, "Blanchard_pm_logFC_scatter_Oligo.3_select.png"), height = 6, width = 10)



