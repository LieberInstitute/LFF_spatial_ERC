## load Blanchard et al. 2022 data from select supp tables

library("tidyverse")
library("here")

#### S12 ####
# Supplementary Table S12. Differentially expressed genes in APOE3/3 vs APOE4/4 iPSC-derived oligodendroglia
list.files(here("external-data", "Blanchard2022"))

Blancard_DE_ipsc_Oligo <- read_csv(here("external-data", "Blanchard2022","Blanchard22_SuppTable12.csv")) |>
    mutate(logFC =  as.numeric(`log2.fold_change.`),
           abs_logFC = abs(logFC))

# Blancard_DE_ipsc_Oligo |> count(q_value < 0.05, status)
# 
# Blancard_DE_ipsc_Oligo |> 
#     filter(status == "OK", is.na(logFC), q_value < 0.05)

Blancard_DE_ipsc_Oligo_signif <- Blancard_DE_ipsc_Oligo |> 
    filter(status == "OK" & q_value < 0.05 & (abs_logFC > 1 | is.na(logFC))) |>
    mutate(reg = ifelse((logFC > 0 | is.na(logFC)), "ipsc_up", "ipsc_down")) |> 
    select(gene, padj = q_value, logFC, reg) 

# Blancard_DE_ipsc_Oligo_signif |> group_by(reg) |> arrange(q_value) |> slice(1:5)
# Blancard_DE_ipsc_Oligo_signif |> count(reg)

# reg           n
# <chr>     <int>
# 1 ipsc_down   662
# 2 ipsc_up    1529

#### S13 ####
# Supplementary Table S13. Differentially expressed genes in post-mortem oligodendrocytes computed by Wilcoxon

grp_table <- tibble(grp = c("APOE34 & APOE44 vs APOE33 (AD and nonAD)", "APOE34 vs APOE33 (AD only)", "APOE34 vs APOE33 (nonAD only)"),
                    grp2 = c("E4+_E33_ALL", "E34_E33_AD", "E34_E33_noAD")) 

Blancard_DE_pm_Oligo <- read_csv(here("external-data", "Blanchard2022","Blanchard22_SuppTable13.csv")) |>
    left_join(grp_table) |>
    dplyr::rename(gene = feature) |>
    mutate(reg = ifelse(logFC > 0, "up", "down"),
           abs_logFC = abs(logFC))

# Blancard_DE_pm_Oligo |> filter(padj < 0.05) |> count(grp2, reg) 

Blancard_DE_pm_Oligo_signif <-  Blancard_DE_pm_Oligo |> 
    filter(padj < 0.05) |> 
    mutate(reg = paste0(grp2, "_", reg)) |> 
    select(gene, padj, logFC, reg) 

# Blancard_DE_pm_Oligo_signif |> count(abs_logFC > 1)
# Blancard_DE_pm_Oligo_signif |> count(reg)

# reg                   n
# <chr>             <int>
# 1 E34_E33_AD_down    1091
# 2 E34_E33_AD_up      2117
# 3 E34_E33_noAD_down  2815
# 4 E34_E33_noAD_up    9467
# 5 E4+_E33_ALL_down   2775
# 6 E4+_E33_ALL_up     8733

#### S14 ####
# Supplementary Table S14. Differential expressed genes in post-mortem oligodendrocytes computed by Nebula

Blancard_DE_Nebula <- read_csv(here("external-data", "Blanchard2022","Blanchard22_SuppTable14.csv")) |>
    dplyr::rename(gene = `...1`, logFC = logFC_Apoe_e4yes) |>
    mutate(reg = ifelse(logFC > 0, "up", "down"),
           abs_logFC = abs(logFC))

# Blancard_DE_Nebula |> filter(padj < 0.05) |> count(reg) 

Blancard_DE_Nebula_signif <- Blancard_DE_Nebula |>
    filter(padj < 0.05) |>
    mutate(reg = paste0("Nebula_", reg)) |> 
    select(gene, padj, logFC, reg) 

# Blancard_DE_Nebula_signif |> count(reg)
# reg             n
# <chr>       <int>
# 1 Nebula_down   528
# 2 Nebula_up     896

#### S15 ####
# Supplementary Table S15. Pseudo-bulk differentially expressed gene statistics for all post-mortem cell types
# 
# Blancard_DE_pb <- read_csv(here("external-data", "Blanchard2022","Blanchard22_SuppTable15.csv")) |>
#     mutate(reg = ifelse(logFC > 0, "up", "down"),
#            abs_logFC = abs(logFC))
# 
# Blancard_DE_pb |> count(adj.P.Val < 0.05, celltype)
# 

#### get gene list ####

Blanchard_DE_All_signif <- Blancard_DE_ipsc_Oligo_signif |>
    bind_rows(Blancard_DE_pm_Oligo_signif) |>
    bind_rows(Blancard_DE_Nebula_signif)

Blanchard_DEG_list <- map(rafalib::splitit(Blanchard_DE_All_signif$reg), ~Blanchard_DE_All_signif$gene[.x])
saveRDS(Blanchard_DEG_list, file = here("external-data", "Blanchard2022", "Blanchard_DEG_list.rds"))

# map_int(Blanchard_DEG_list, length)
# E34_E33_AD_down     E34_E33_AD_up E34_E33_noAD_down   E34_E33_noAD_up  E4+_E33_ALL_down    E4+_E33_ALL_up         ipsc_down           ipsc_up       Nebula_down 
# 1091              2117              2815              9467              2775              8733               662              1529               528 
# Nebula_up 
# 896 
