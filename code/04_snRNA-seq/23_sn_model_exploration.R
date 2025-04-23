## April 2025, Louise Huuki-Myers
## Explore modeling results from sn modeling

library("spatialLIBD")
library("purrr")
library("tidyverse")
library("ComplexHeatmap")
library("ggrepel")
library("here")
library("sessioninfo")

#### Set up dirs ####
# data_dir <- here("processed-data", "04_snRNA-seq", "23_sn_model_exploration")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "23_sn_model_exploration")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load sn ERC modeling ####
erc_sn_modeling_results <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk", "modeling_results-cell_type_anno.rds"))


colnames(erc_sn_modeling_results$pairwise)
# Oligo.1-Oligo.2
# Astro.1-Astro.2
# Micro.1-Micro.2

##### Microglia ####

micro_pairwise <- erc_sn_modeling_results$pairwise |> 
    select(ensembl, gene, ends_with("Micro.1-Micro.2")) |>
    rename_with(~gsub("_Micro.1-Micro.2", "", .x))  |>
    mutate(DE_class = case_when(logFC > 1 & fdr < 0.05 ~ "Micro.1",
                                logFC < -1 & fdr < 0.05 ~ "Micro.2",
                                TRUE ~"Other"))


micro_pairwise |>
    filter(DE_class != "Other") |>
    group_by(DE_class) |>
    # arrange(-abs(logFC)) |>
    arrange(fdr) |>
    slice(1:5)

# ensembl         gene            t_stat  p_value      fdr logFC DE_class
# <chr>           <chr>            <dbl>    <dbl>    <dbl> <dbl> <chr>   
# 1 ENSG00000100351 GRAP2            -5.32 1.53e- 7 7.07e- 7 -1.72 Micro.1 
# 2 ENSG00000127152 BCL11B           -4.91 1.23e- 6 4.55e- 6 -1.41 Micro.1 
# 3 ENSG00000141293 SKAP1            -4.65 4.25e- 6 1.39e- 5 -1.83 Micro.1 
# 4 ENSG00000227191 TRGC2            -4.29 2.18e- 5 6.04e- 5 -1.07 Micro.1 
# 5 ENSG00000172673 THEMIS           -3.78 1.75e- 4 4.01e- 4 -1.48 Micro.1 
# 6 ENSG00000275799 ENSG00000275799  15.7  1.74e-45 4.51e-41  2.77 Micro.2 
# 7 ENSG00000188343 CIBAR1           14.0  1.09e-37 1.41e-33  2.19 Micro.2 
# 8 ENSG00000288900 ENSG00000288900  13.8  9.70e-37 8.39e-33  2.06 Micro.2 
# 9 ENSG00000196866 H2AC7            13.7  2.20e-36 1.42e-32  2.16 Micro.2 
# 10 ENSG00000165383 LRRC18           13.1  7.06e-34 3.67e-30  1.43 Micro.2 


micro_pairwise |> filter(gene %in% c("CST7", "LPL"))

# Sestan clusters
# Micro.C1QB.P2RY12 - Micro2
# Micro.C1QB.CD83 - Micro1 ?
micro_pairwise |> filter(gene %in% c("P2RY12", "CD83"))
#                         ensembl   gene   t_stat      p_value         fdr     logFC DE_class
# ENSG00000169313 ENSG00000169313 P2RY12 1.423956 0.1550824673 0.196568050 0.4035809    Other
# ENSG00000112149 ENSG00000112149   CD83 3.426462 0.0006618343 0.001346388 0.9442603    Other


(pval_lim <- max(micro_pairwise$p_value[micro_pairwise$fdr < 0.05]))

micro_colors <- c(cell_type_colors$anno[grepl("Micro", names(cell_type_colors$anno))],
                  other = "grey")

pairwise_volcano_Micro <- micro_pairwise |>
    ggplot(aes(`logFC`, -log10(`p_value`), color = DE_class)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = ifelse(DE_class != "Other", gene, "")))  +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed") +
    scale_color_manual(values = micro_colors)

ggsave(pairwise_volcano_Micro, filename = here(plot_dir, "pairwise_volcano_Micro.png"))

##### Oligo ####
oligo_pairwise <- erc_sn_modeling_results$pairwise |> 
    select(ensembl, gene, ends_with("Oligo.1-Oligo.2")) |>
    rename_with(~gsub("_Oligo.1-Oligo.2", "", .x))  |>
    mutate(DE_class = case_when(logFC > 1 & fdr < 0.05 ~ "Oligo.1",
                                logFC < -1 & fdr < 0.05 ~ "Oligo.2",
                                TRUE ~"Other"))

# Sestan clusters
# Oligo.OPALIN.SLC5A11 (Oligo.2)
# Oligo.OPALIN.LAMA2 (Oligo.1)
# Oligo.OPALIN.LINC01098
oligo_pairwise |> filter(gene %in% c("SLC5A11", "LAMA2","LINC01098"))
#                      ensembl      gene      t_stat    p_value         fdr       logFC DE_class 
# ENSG00000231171 ENSG00000231171 LINC01098 -0.866699 0.386522272 0.59134928 -0.2730363    Other X
# ENSG00000196569 ENSG00000196569     LAMA2  2.536581 0.011496825 0.04799157  0.7007831    Other Oligo.1
# ENSG00000158865 ENSG00000158865   SLC5A11 -3.135256 0.001817941 0.01050222 -0.9262855    Other Oligo.2

oligo_colors <- c(cell_type_colors$anno[grepl("Oligo", names(cell_type_colors$anno))],
                  other = "grey")

(pval_lim <- max(oligo_pairwise$p_value[oligo_pairwise$fdr < 0.05]))

pairwise_volcano_Oligo <- oligo_pairwise |>
    ggplot(aes(`logFC`, -log10(`p_value`), color = DE_class)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = ifelse(DE_class != "Other", gene, "")))  +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed") +
    scale_color_manual(values = oligo_colors)

ggsave(pairwise_volcano_Oligo, filename = here(plot_dir, "pairwise_volcano_Oligo.png"))

##### Astro ####
astro_pairwise <- erc_sn_modeling_results$pairwise |> 
    select(ensembl, gene, ends_with("Astro.1-Astro.2")) |>
    rename_with(~gsub("_Astro.1-Astro.2", "", .x))  |>
    mutate(DE_class = case_when(logFC > 1 & fdr < 0.05 ~ "Astro.1",
                                logFC < -1 & fdr < 0.05 ~ "Astro.2",
                                TRUE ~"Other"))

# Sestan clusters
# Astro.AQP4.GFAP
# Astro.AQP4.CHRDL1

astro_pairwise |> filter(gene %in% c("GFAP", "CHRDL1"))
#                    ensembl      gene      t_stat    p_value         fdr       logFC DE_class
# ENSG00000131095 ENSG00000131095   GFAP -2.465419 0.01401978 0.0461272 -0.7459056    Other (Astro.2)
# ENSG00000101938 ENSG00000101938 CHRDL1  1.315948 0.18879404 0.3243615  0.3591450    Other (Astro.1)

astro_colors <- c(cell_type_colors$anno[grepl("Astro", names(cell_type_colors$anno))],
                  other = "grey")

(pval_lim <- max(astro_pairwise$p_value[astro_pairwise$fdr < 0.05]))

pairwise_volcano_Astro <- astro_pairwise |>
    ggplot(aes(`logFC`, -log10(`p_value`), color = DE_class)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = ifelse(DE_class != "Other", gene, "")))  +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed") +
    scale_color_manual(values = astro_colors)

ggsave(pairwise_volcano_Astro, filename = here(plot_dir, "pairwise_volcano_Astro.png"))
