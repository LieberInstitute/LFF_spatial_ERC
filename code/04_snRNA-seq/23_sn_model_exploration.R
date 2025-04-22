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


##### Microglia ####

colnames(erc_sn_modeling_results$pairwise)

micro_pairwise <- erc_sn_modeling_results$pairwise |> 
    select(ensembl, gene, ends_with("Micro.1-Micro.2")) |>
    rename_with(~gsub("_Micro.1-Micro.2", "", .x))  |>
    mutate(DE_class = case_when(logFC > 1 & fdr < 0.05 ~ "Micro.2",
                                logFC < -1 & fdr < 0.05 ~ "Micro.1",
                                TRUE ~"Other"))


micro_pairwise |>
    filter(DE_class != "Other") |>
    group_by(DE_class) |>
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

(pval_lim <- max(micro_pairwise$p_value[micro_pairwise$fdr < 0.05]))

micro_colors <- c(cell_type_colors$anno[grepl("Micro", names(cell_type_colors$anno))],
                  other = "grey")

micro_pairwise |>
    ggplot(aes(`logFC`, -log10(`p_value`), color = DE_class)) +
    geom_point(size = 1) +
    geom_text_repel(aes(label = ifelse(DE_class != "Other", gene, "")))  +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_lim), linetype = "dashed") +
    scale_color_manual(values = micro_colors)

