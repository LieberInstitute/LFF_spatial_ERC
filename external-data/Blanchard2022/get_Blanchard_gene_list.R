Blanchard_gene_list <- list(Blanchard_Fig2c = c("DHCR24", "LPCAT3", "SCP2", "PRKN","SEC23A","PCYT1B","MBTPS1","LPIN1","LBR","NR1H2","IRS2","LPIN2"),
                           Blanchard_Fig4a = c("PLLP", "MYRF", "MAG", "OPALIN", "SREBF1", "PLP1", "SCAP", "MVK", "FDPS", "ABCG1", "HMGCS1", "IDI1", "LDLR", "INSIG1", "SREBF2", "SQLE", "DHCR7", "DHCR24", "FDFT1", "LSS"),
                           Blanchard_ExFig8a = c("PLP1", "OPALIN", "PLLP", "MYRF", "MAG", "MOG"))

blanchard_gene_annotations <- map(Blanchard_gene_list, ~tibble(gene_name = .x, anno = ""))

walk2(blanchard_gene_annotations, names(blanchard_gene_annotations), ~write_csv(.x, paste0(.y, "_gene_anno.csv")))


## combine w/ Wilcoxon test stats
Blanchard_gene_annotaions <- map_dfr(list.files(pattern = "Fig"), read_csv) |> unique()

# Supplementary Table S13. Differentially expressed genes in post-mortem oligodendrocytes computed by Wilcoxon

grp_table <- tibble(grp = c("APOE34 & APOE44 vs APOE33 (AD and nonAD)", "APOE34 vs APOE33 (AD only)", "APOE34 vs APOE33 (nonAD only)"),
                    grp2 = c("E4+_E33_ALL", "E34_E33_AD", "E34_E33_noAD")) 

Blancard_DE_pm_Oligo <- read_csv(here("external-data", "Blanchard2022","Blanchard22_SuppTable13.csv")) |>
    left_join(grp_table) |>
    dplyr::rename(gene = feature) |>
    mutate(reg = ifelse(logFC > 0, "up", "down"),
           abs_logFC = abs(logFC))

Blanchard_gene_anno_stats <- Blanchard_gene_annotaions |>
    left_join(Blancard_DE_pm_Oligo |> select(gene_name = gene, grp2, logFC, padj)) |>
    mutate(p.signif = case_when(padj < 0.005 ~ "***",
                                           padj < 0.01 ~"**",
                                           padj < 0.05 ~"*",
                                           TRUE~""))

Blanchard_gene_anno_stats |> count(p.signif)

Blanchard_stats_col <- Blanchard_gene_anno_stats |>
    ggplot(aes(y = gene_name, x = logFC, fill = anno)) +
    geom_col() +
    facet_grid(anno~grp2, scales = "free", space = "free") +
    # geom_text(aes(label = p.signif)) +
    theme_bw()

ggsave(Blanchard_stats_col, filename = "Blanchard_stats_col.png")
