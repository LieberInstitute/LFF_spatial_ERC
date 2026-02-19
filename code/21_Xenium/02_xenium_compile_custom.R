## Louise Huuki-Myers, Jan 2026
## Build custom list

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("readxl")
library("writexl")
library("ggrepel")

data_dir <- here("processed-data", "21_Xenium", "02_xenium_compile_custom")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "02_xenium_compile_custom")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load and annotate Xenium base panel ####
Xenium_annotation <- read_xlsx(here("processed-data", "21_Xenium", "01_xenium_panel", "Xenium_hBrain_annotation.xlsx"))
Xenium_annotation |> dplyr::count(cell_type_class)

Xenium_hBrain <- read_csv(here("processed-data", "21_Xenium", "Xenium_hBrain_v1_metadata.csv")) |>
    left_join(Xenium_annotation) |>
    arrange(anno) 

#### Lit - AD risk genes ####
AD_risk <- read_csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv"))

Xenium_hBrain |> dplyr::filter(Genes %in% AD_risk$symbol) |> select(Genes, Annotation, cell_type_anno)
# Genes  Annotation                  cell_type_anno
# <chr>  <chr>                       <chr>         
# 1 APP    Broad                       NA            
# 2 SORCS1 L5 ET                       Excit.L5      
# 3 PSEN2  Glioblastoma (Cancer cells) NA            
# 4 APOE   Microglia-PVM               NA            
# 5 PSEN1  Oligodendrocyte             NA  

## select AD risk genes already in panel 
select_AD_risk <- Xenium_hBrain |> 
    dplyr::filter(Genes %in% AD_risk$symbol) |> 
    transmute(gene_name = Genes, 
              target = NA, 
              goal = "Lit - AD risk", 
              note = "OpenTargets AD gene", 
              in_base = TRUE)

#### Lit - NE_receptor_genes ####

NE_receptor_genes <- c("ADRA1A","ADRA1B","ADRA1D","ADRA2A","ADRA2B","ADRA2C","ADRB1","ADRB2","ADRB3")

# pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_NE_receptor.pdf")), height = 7, width = 8)
# sce |>
#     scDotPlot(features = NE_receptor_genes,
#               group = "cell_type_anno",
#               groupAnno = "cell_type_anno",
#               scale = TRUE,
#               annoColors = list("cell_type_anno" = cell_type_colors$anno),
#               clusterRows = FALSE,
#               groupLegends = FALSE
#     ) 
# dev.off()

## select genes expressed in Oligo & Astro
select_NE_receptors <- tibble(gene_name = c("ADRA1A","ADRA1B","ADRB1"), 
                              target = "Oligo, Astro", 
                              goal = "Lit - NE receptors") |>
    mutate(in_base = gene_name %in% Xenium_hBrain$Genes)


#### Marker Genes - Cell Type  - global enrichment ####

## global enrichment didn't look great, did not include

enrichment_stats_top <- read_csv(here("processed-data", "04_snRNA-seq", "31_sn_subcluster_heatmap", "sce_subcluster_enrichment_top5.csv")) |>
    dplyr::rename(Genes = ensembl)

non_unique <- enrichment_stats_top |> count(Genes) |> filter(n != 1) |> pull(Genes)

# enrichment_stats_top |> inner_join(Xenium_hBrain)

any(enrichment_stats_top |> filter(Genes %in% Xenium_hBrain$Genes)) # no overlaps 

enrichment_stats_select <- enrichment_stats_top |> 
    filter(cell_type_anno %in% c("Astro.1", "Astro.2", "Astro.3",
                                 "Oligo.1", "Oligo.2", "Oligo.3", "Oligo.5",
                                 "OPC.5", 
                                 "Excit.L5.2", "Excit.L2_5.1")
    ) |> 
    dplyr::rename(gene_name = Genes)|>
    filter(gene_name != "ENSG00000272076") ## top gene for several astros


# rowData(sce)$enrichment <- NULL
# rowData(sce)$enrichment <- enrichment_stats_select$cell_type_anno[match(rownames(sce), enrichment_stats_select$gene_name)] 
# table(rowData(sce)$enrichment)
# 
# pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_enrichment_select.pdf")), height = 11, width = 8)
# sce |>
#     scDotPlot(features = enrichment_stats_select$gene_name,
#               group = "cell_type_anno",
#               groupAnno = "cell_type_anno",
#               featureAnno = "enrichment",
#               scale = TRUE,
#               annoColors = list("cell_type_anno" = cell_type_colors$anno,
#                                 "enrichment" = cell_type_colors$anno),
#               clusterRows = FALSE,
#               groupLegends = FALSE
#     ) 
# dev.off()

#### Marker Genes - Cell Type - Global MeanRatio ####

load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio","marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats_MeanRatio |>
    filter(MeanRatio.rank <= 10, MeanRatio > 1) |>
    arrange(cellType.target)

## 21 top 10 genes overlap
mean_ratio_xenium_overlap <- marker_stats_top |> 
    inner_join(Xenium_hBrain |> 
                   dplyr::rename(gene = Genes)) |>  
    select(gene, cellType.target, MeanRatio, MeanRatio.rank, Annotation, cell_type_anno, MeanRatio.anno) |>
    arrange(cellType.target) 

mean_ratio_xenium_overlap |> print(n = 21)

## pick markers for cell types of interest - not covered by base

cellTypes_of_interest <- c("Astro.1", "Astro.2", "Astro.3",
                           "Oligo.1", "Oligo.2", "Oligo.3", "Oligo.5",
                           "OPC.5",
                           "Excit.L5.2", "Excit.L2_5.1")

marker_stats_select <- marker_stats_top |> 
    filter(MeanRatio.rank <= 5) |>
    filter(cellType.target %in% cellTypes_of_interest
    )

## filter to very specific genes
marker_stats_select <- marker_stats_select |> filter(MeanRatio > 1.5, MeanRatio.rank <= 2)

## select target global MR genes or MR genes overlapping base panel
select_global_mean_ratio_marker <- mean_ratio_xenium_overlap |>
    bind_rows(marker_stats_select) |>
    transmute(gene_name = gene, 
              target = cellType.target, 
              goal = "Marker - cell type global",
              note = paste("globl MR genes - ", MeanRatio.anno),
              in_base = gene_name %in% Xenium_hBrain$Genes) |>
    arrange(target)

select_global_mean_ratio_marker |> dplyr::count(target, in_base)
# target           in_base     n
# <fct>            <lgl>   <int>
# 1 Macro            TRUE        3
# 2 Micro.1          TRUE        1
# 3 OPC.4            TRUE        1
# 4 Oligo.1          TRUE        1
# 5 Oligo.2          FALSE       1
# 6 Vasc.VLMC        TRUE        1
# 7 Excit.L5.2       FALSE       1
# 8 Excit.L5_6_NP    TRUE        1
# 9 Excit.L6_CT      TRUE        1
# 10 Excit.L6b        TRUE        2
# 11 Inhib.Chandelier TRUE        1
# 12 Inhib.Lamp5_Lhx6 TRUE        3
# 13 Inhib.Pvalb      TRUE        2
# 14 Inhib.Sst        TRUE        2
# 15 Inhib.Vip        TRUE        2


# rowData(sce)$MeanRatio <- NULL
# rowData(sce)$MeanRatio <- as.character(select_global_mean_ratio_marker$target)[match(rownames(sce), select_global_mean_ratio_marker$gene_name)] 
# table(rowData(sce)$MeanRatio)
# 
# pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_MeanRatio_global_select.pdf")), height = 11, width = 8)
# sce |>
#     scDotPlot(features = select_global_mean_ratio_marker$gene_name,
#               group = "cell_type_anno",
#               groupAnno = "cell_type_anno",
#               featureAnno = "MeanRatio",
#               scale = TRUE,
#               annoColors = list("cell_type_anno" = cell_type_colors$anno,
#                                 "MeanRatio" = cell_type_colors$anno),
#               clusterRows = FALSE,
#               clusterColumns = FALSE,
#               groupLegends = FALSE
#     ) 
# sce |>
#     scDotPlot(features = select_global_mean_ratio_marker$gene_name,
#               group = "cell_type_anno",
#               groupAnno = "cell_type_anno",
#               featureAnno = "MeanRatio",
#               scale = TRUE,
#               annoColors = list("cell_type_anno" = cell_type_colors$anno,
#                                 "MeanRatio" = cell_type_colors$anno),
#               clusterRows = TRUE,
#               clusterColumns = TRUE,
#               groupLegends = FALSE
#     ) 
# dev.off()

#### Marker Genes - Cell Type - cell type specific modeling ####

broad_cell_types <- c("Oligo", "Astro", "OPC", "Excit")
names(broad_cell_types) <- broad_cell_types

ct_specific_mod <- map(broad_cell_types, ~readRDS(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", .x, sprintf("modeling_results_subtype-%s.rds", .x))) |>
                           pluck("enrichment") |>
                           pivot_longer(!c("ensembl", "gene"), names_to = "stat") |>
                           mutate(stat = gsub("p_value", "p.value", stat),
                                  stat = gsub("t_stat", "t.stat", stat),) |>
                           separate(stat, into = c("stat", "cellType.target"), sep = "_", extra = "merge") |>
                           pivot_wider(names_from = "stat", values_from = "value", names_prefix = "enrich_") |>
                           dplyr::rename(gene_ensembl = ensembl)
)

head(ct_specific_mod$Oligo)

ct_specific_mod$Oligo |> filter(gene == "LINGO2")

ct_specific_mod$Oligo |> filter(cellType.target == "Oligo.3") |> arrange(-enrich_t.stat)

ct_specific_marker_stats <- map2(broad_cell_types, ct_specific_mod, 
                                 ~get(load(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", .x, sprintf("marker_stats_MeanRatio_%s.Rdata", .x)))) |>
                                     full_join(.y) |>
                                     mutate(in_base = gene %in% Xenium_hBrain$Genes))

ct_specific_marker_stats$Oligo |> filter(gene == "GPM6A")

mean_ratio_specific_xenium_overlap <- map_dfr(ct_specific_marker_stats, 
                                              ~.x |> 
                                                  filter(MeanRatio > 1.5 & in_base) |>
                                                  transmute(gene_name = gene, 
                                                            target = cellType.target, 
                                                            goal = "Marker - cell type specific",
                                                            note = paste("specific MR genes -", MeanRatio.anno),
                                                            in_base = in_base))



enrichment_specific_xenium_overlap <- map_dfr(ct_specific_marker_stats,
                                          ~.x |>
                                              filter(enrich_fdr < 0.05)|>
                                              group_by(cellType.target) |>
                                              arrange(-enrich_t.stat) |>
                                              slice(1:10) |>
                                              filter(in_base) |>
                                              transmute(gene_name = gene,
                                                        target = cellType.target,
                                                        goal = "Marker - cell type specific",
                                                        note = paste("specific Enrich genes t=", round(enrich_t.stat, 2)),
                                                        in_base = in_base))

enrichment_specific_xenium_overlap |> count(target)


walk2(ct_specific_marker_stats, names(ct_specific_marker_stats), function(marker_stats, name){
    
    ct_specific_t_v_ratio <- marker_stats |>
        ggplot(aes(x= MeanRatio, y = enrich_t.stat, color = enrich_fdr < 0.05)) +
        geom_point(size = 0.5) +
        scale_color_manual(values = c(`TRUE` = "red")) +
        geom_text_repel(aes(label = gene), size = 1.2) +
        facet_wrap(~cellType.target, scales = "free_x") +
        theme_bw() +
        geom_vline(xintercept = 1, color = "blue") +
        theme(legend.position = "None")
    
    ggsave(ct_specific_t_v_ratio, filename = here(plot_dir, sprintf("ct_specific_t_v_ratio_%s.png", name)))
    
})

mean_ratio_specific_interest <- map_dfr(ct_specific_marker_stats, 
                                        ~.x |> 
                                            filter(cellType.target %in% cellTypes_of_interest,
                                                   MeanRatio > 1, 
                                                   MeanRatio.rank <= 3) |>
                                            mutate(note = paste0("specific MR genes - ", MeanRatio.anno)))

mean_ratio_specific_interest  |> dplyr::count(cellType.target)

enrichment_specific_interest <- map_dfr(ct_specific_marker_stats, 
                                        ~.x |> 
                                            filter(cellType.target %in% cellTypes_of_interest,
                                                   enrich_p.value < 0.05, 
                                                   enrich_t.stat > 0) |>
                                            arrange(-enrich_t.stat) |>
                                            group_by(cellType.target) |>
                                            dplyr::slice(1:5) |>
                                            mutate(note = paste0("specific Enrich genes t=", round(enrich_t.stat, 2))))

ct_specific_marker_stats$Oligo |>
    filter(cellType.target == "Oligo.3", gene %in% c("LINGO2"))


enrichment_specific_interest |> filter(cellType.target == "Oligo.3")

enrichment_specific_interest|> dplyr::count(cellType.target)

## TODO edit down n markers for non Oligo.3 cell types
select_specific_cell_type_marker <- bind_rows(mean_ratio_specific_interest, 
                                              enrichment_specific_interest)  |>
    transmute(gene_name = gene, 
              target = cellType.target, 
              goal = "Marker - cell type specific",
              note = note,
              in_base = in_base) |>
    bind_rows(mean_ratio_specific_xenium_overlap) |>
    bind_rows(enrichment_specific_xenium_overlap) |>
    group_by(gene_name, target, goal, in_base) |>
    summarise(note = paste(unique(note), collapse = ", "))

select_specific_cell_type_marker |> ungroup() |> filter(!in_base) |> dplyr::count(target)

# target           n
# <chr>        <int>
#     1 Astro.1          5
# 2 Astro.2          4
# 3 Astro.3          5
# 4 Excit.L2_5.1     5
# 5 Excit.L5.2       5
# 6 OPC.5            4
# 7 Oligo.1          4
# 8 Oligo.2          4
# 9 Oligo.3          4
# 10 Oligo.5          2

select_specific_cell_type_marker |> ungroup() |> count(target)

select_specific_cell_type_marker |> filter(target == "Oligo.3")
select_specific_cell_type_marker |> filter(target == "Astro.1")
select_specific_cell_type_marker |> filter(target == "Astro.2")

ct_specific_marker_stats$Oligo |> filter(cellType.target == "Oligo.3", gene_name %in% c("LINGO2", "KCNJ3", "GPM6A", "MT-CO3", "CNTNAP2"))

# ct_markers <- c(
#     Oligo.1_2 <- c("OPALIN"),
#     Oligo.3 <- c("LINGO2", "KCNJ3", "GPM6A", "MT-CO3", "CNTNAP2"),
#     Oligo.5 <- c("ARHGEF3")
# )

#### Marker Gene - SpD ####

load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_SpD.Rdata"), verbose = TRUE) 

SpD_marker_stats <- marker_stats |>
    mutate(in_base =  gene %in% Xenium_hBrain$Genes)

# SpD_marker_stats |> filter(cellType.target == "LD~Sp09D02")

## Get top enrichment genes for WM.uf 
SpD_enrichment <- read_csv(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "enrichment_modeling_SpD_top100.csv")) |>
    mutate(in_base =  gene %in% Xenium_hBrain$Genes)  |> 
    filter((top <= 5 & in_base) | (top <= 2 & test == "WM.uf~Sp09D07")) |>
    transmute(gene_name = gene, 
              target = test, 
              goal = "Marker - SpD",
              note = paste0("SpD Enrich gene t=", round(stat, 2)),
              in_base = in_base)

## top MR genes for LD
select_SpD_markers <- SpD_marker_stats |> 
    filter(MeanRatio >1 & ((MeanRatio.rank <= 5 & in_base) | (MeanRatio.rank <= 2 & cellType.target == "LD~Sp09D02"))) |>
    transmute(gene_name = gene, 
              target = cellType.target, 
              goal = "Marker - SpD",
              note = paste0("SpD MeanRatio gene: ", MeanRatio.anno),
              in_base = in_base) |>
    bind_rows(SpD_enrichment) |>
    group_by(gene_name, target, goal, in_base) |>
    summarise(note = paste0(note, collapse = ", ")) |>
    ungroup() |>
    add_row(gene_name = "RELN", ## Add RELN as ERC Layer 2 marker
            goal = "Marker - SpD",
            target = "L2.3~Sp09D01",
            note = paste0("SpD marker from Lit"),
            in_base = "RELN" %in% Xenium_hBrain$Genes)

select_SpD_markers |> dplyr::count(gene_name, target) |> arrange(-n)


#### DEGs - Oligo ####

#### DEGs - Excit ####

#### DEGs - Astro Mediation Genes ####

# list.files(here("processed-data","22_Mediation","out-erc_astro"))

mediation_summary <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediation_vs_baseline_summary.tsv.gz"))

mediator_tab <- mediation_summary |> 
    filter(mediated_n > 0) |> 
    select(mediation_run, mediated_n) |> 
    separate(mediation_run, into = c("cell_type", "ensemblID", "gene_name"), sep = "\\|") |>
    mutate(in_base = gene_name %in% Xenium_hBrain$Genes)
# cell_type       ensemblID gene_name mediated_n
# 1   Astro.1 ENSG00000117600    PLPPR4          1
# 2   Astro.1 ENSG00000221890     NPTXR         22
# 3   Astro.2 ENSG00000141338     ABCA8          1
# 4   Astro.2 ENSG00000147488      ST18          1
# 5   Astro.2 ENSG00000166342     NETO1         61
# 6   Astro.2 ENSG00000177283      FZD8        125
# 7   Astro.2 ENSG00000185518      SV2B        199
# 8   Astro.3 ENSG00000133019     CHRM3        205
# 9   Astro.3 ENSG00000134352     IL6ST          4

mediated_gene_shortlist <- c("MBP", # mylenation - ABCA8 outcome
                             "NAP1L3", # only outcome pf ST18
                             "ENC1", # only outcome of PLPPR4
                             "SOX6", # key regulator of oligodendrocyte maturation - NPTXR outcome
                             "SORBS1", # IL6ST outcome
                             "GPM6A" #Oligo.3 marker FZD8 outcome
)

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(outcome_in_base = outcome %in% Xenium_hBrain$Genes,
           outcome_in_list = outcome %in% mediated_gene_shortlist,
           pair = paste0(mediator, "|", outcome))

## 15 pairs from short list
mediator_outcome |> filter(outcome_in_list) 

## 20 outcomes in base probes
mediator_outcome |> filter(outcome_in_base) 

mediator_outcome |> filter(outcome_in_list | outcome_in_base) 

mediator_outcome |> filter(outcome_in_base) |> dplyr::count(med_cl, mediator)

mediator_outcome |> dplyr::count(med_cl, mediator)

mediator_outcome |> filter(mediator == "FZD8")

mediation_fdr_scatter <- mediator_outcome |> 
    ggplot(aes(x = fdr, fdr_med, color = outcome_in_base)) +
    geom_point() +
    geom_text_repel(aes(label = outcome), size = 1.5) +
    facet_wrap(~mediator) +
    theme_bw()

ggsave(mediation_fdr_scatter, filename = here(plot_dir, "mediation_fdr_scatter.png"), height = 8, width = 10)

mediation_logFC_scatter <- mediator_outcome |> 
    ggplot(aes(x = logFC, logFC_med, color = outcome_in_base)) +
    geom_point() +
    geom_text_repel(aes(label = outcome), size = 1.5) +
    facet_wrap(~mediator) +
    theme_bw()

ggsave(mediation_logFC_scatter, filename = here(plot_dir, "mediation_logFC_scatter.png"), height = 8, width = 10)

## select all mediator genes, and outcome genes in shortlist or in base panel, 
mediator_selection <- mediator_outcome |> 
    select(mediator, outcome, med_cl) |> 
    pivot_longer(!c(med_cl), values_to = "gene_name", names_to = "role") |>
    mutate(in_base = gene_name %in% Xenium_hBrain$Genes,
           in_list = gene_name %in% mediated_gene_shortlist) |>
    filter(in_base | in_list | grepl("mediator", role)) |>
    group_by(gene_name, in_base, in_list) |>
    summarise(med_cl = paste(unique(med_cl),collapse = ","),
              role = paste(unique(role),collapse = ",")) 

mediator_selection |> ungroup() |> dplyr::count(in_base)

## select mediation genes
select_mediator_genes <- mediator_selection |>
    mutate(target = ifelse(role == "outcome", "Oligo.3", med_cl),
           goal = "DEG - mediation", 
           note = paste("mediation gene:", role, target)) |>
    ungroup() |>
    select(gene_name, in_base,target, goal, note)

mediator_tested <- mediator_outcome |> 
    filter(mediator %in% mediator_selection$gene_name & outcome %in% mediator_selection$gene_name)

mediator_tested |> dplyr::count(med_cl)

mediator_selection |> filter(role == "outcome") |> arrange(-n_pairs)

mediator_selection |>
    dplyr::count(in_base, in_list, role)


pdf(here(plot_dir, "sn_dotplot_mediator_genes.pdf"))
sce |>
    scDotPlot(features = mediator_selection_add$gene_name,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = FALSE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 

sce |>
    scDotPlot(features = mediator_selection_add$gene_name,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 
dev.off()

#### COMPILE ALL PROBE LISTS ####

ALL_probes_long <- select_AD_risk |>
    bind_rows(select_global_mean_ratio_marker) |>
    bind_rows(select_specific_cell_type_marker) |>
    bind_rows(select_SpD_markers) |>
    bind_rows(select_NE_receptors) |>
    bind_rows(select_mediator_genes)

# ALL_probes_long |> filter(goal == "")

ALL_probes_summary <- ALL_probes_long |>
    group_by(gene_name, in_base) |>
    summarise(n_goal = n(), 
              goals = paste(sort(unique(goal)), collapse = ", "))

# ALL_probes_summary |> filter(goals == "")

probes_summary_count <- ALL_probes_summary |>
    ungroup() |>
    dplyr::count(goals, in_base) |>
    pivot_wider(names_from = "in_base", names_prefix = "base_", values_from = "n") 


bind_rows(probes_summary_count,
          summarise(probes_summary_count,
                    across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), # Sum all numeric columns
                    across(where(is.character), ~"Total") # Label the character column as "Total"
          )
)

writexl::write_xlsx(ALL_probes_long, path = here(data_dir, "ERC_Xenium_ALL_probes_long.xlsx"))

## Oligo.3 check 

ALL_probes_long |> filter(target == "Oligo.3")


## filter gene not in base
select_custom_probes <- ALL_probes_gene_summary |> filter(!in_base)


