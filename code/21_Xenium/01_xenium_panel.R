## Louise Huuki-Myers, Jan 2026
## Explore existing panel

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("scDotPlot")
library("readxl")
library("ggrepel")

data_dir <- here("processed-data", "21_Xenium", "01_xenium_panel")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "01_xenium_panel")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load sce data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors

## spatial data
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

rownames(spe) <- rowData(spe)$gene_name

SpD_colors <- metadata(spe)$SpD_colors

#### load and annotate Xenium base panel ####
Xenium_annotation <- read_xlsx(here(data_dir, "Xenium_hBrain_annotation.xlsx"))
Xenium_annotation |> dplyr::count(cell_type_class)

Xenium_hBrain <- read_csv(here("processed-data", "21_Xenium", "Xenium_hBrain_v1_metadata.csv")) |>
    mutate(in_sce = Genes %in% rownames(sce),
           in_spe = Genes %in% rownames(spe)) |>
    left_join(Xenium_annotation) |>
    arrange(anno)  |> 
    filter(in_sce & in_spe)

Xenium_hBrain |> filter(!in_sce | !in_spe)
# Genes Ensembl_ID Num_Probesets Codewords Annotation in_sce in_spe anno  cell_type_anno cell_type_broad cell_type_class
# <chr> <chr>              <dbl>     <dbl> <chr>      <lgl>  <lgl>  <chr> <chr>          <chr>           <chr>          
#     1 BTBD… ENSG00000…             8         1 Pvalb      FALSE  TRUE   Inhi… Inhib.Pvalb    Inhib           Neuron  

# Xenium_hBrain <- Xenium_hBrain |> filter(in_sce)

# Xenium_hBrain |> dplyr::count(Annotation) |> write_csv(here(data_dir, "Xenium_hBrain_summary.csv"))

Xenium_hBrain |> dplyr::count(Annotation) |> print(n = 28) 

# Annotation                      n
# <chr>                       <int>
# 1 Astrocyte                      11
# 2 Broad                           1
# 3 Chandelier                      9
# 4 Endothelial                     8
# 5 Glioblastoma (Cancer cells)    24
# 6 Glioblastoma (TME)             44
# 7 L2/3 IT                         3
# 8 L4 IT                           6
# 9 L5 ET                           5
# 10 L5 IT                           2
# 11 L5/6 NP                         7
# 12 L6 CT                           4
# 13 L6 IT                           3
# 14 L6 IT Car3                     11
# 15 L6b                             3
# 16 Lamp5                           2
# 17 Lamp5 Lhx6                      8
# 18 Microglia-PVM                  19
# 19 OPC                            12
# 20 Oligodendrocyte                23
# 21 Pax6                            9
# 22 Proliferation                   8
# 23 Pvalb                           6
# 24 Sncg                            5
# 25 Sst                             4
# 26 Sst Chodl                       8
# 27 VLMC                           18
# 28 Vip                             2

#### AD risk genes ####
AD_risk <- read_csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv"))

Xenium_hBrain |> dplyr::filter(Genes %in% AD_risk$symbol) |> select(Genes, Annotation, cell_type_anno)
# Genes  Annotation                  cell_type_anno
# <chr>  <chr>                       <chr>         
# 1 APP    Broad                       NA            
# 2 SORCS1 L5 ET                       Excit.L5      
# 3 PSEN2  Glioblastoma (Cancer cells) NA            
# 4 APOE   Microglia-PVM               NA            
# 5 PSEN1  Oligodendrocyte             NA  

#### dotplot of panel on cell_type_anno ####

rowData(sce)$Xenium <- NULL
rowData(sce)$Xenium <- Xenium_hBrain$anno[match(rownames(sce), Xenium_hBrain$Genes)] 
table(rowData(sce)$Xenium)

xenium_colors <- c(cell_type_colors$anno[Xenium_annotation$anno[Xenium_annotation$anno %in% names(cell_type_colors$anno)]],
                   cell_type_colors$broad[Xenium_annotation$anno[Xenium_annotation$anno %in% names(cell_type_colors$broad)]],
                   Inhib.Lamp5 = "#9A1D1D",
                   Inhib.Sncg = "#D39E9E",
                   Inhib.Sst_Chold = "#660000",
                   Excit.L2_3_IT = "#7290E9",
                   Excit.L4_IT = "#7DCEC0",
                   Excit.L5_ET = "#70B8FF",
                   Excit.L5_IT = "#295D91",
                   Excit.L6_IT = "#293F91",
                   Excit.L6_IT_Car3 = "#2F7085",
                   Broad = "#BEB7A4",
                   Glioblastoma_cancer = "#705748",
                   Glioblastoma_TME = "#B18859",
                   Proliferation = "#F7FFE0")

Xenium_annotation$anno[!Xenium_annotation$anno %in% names(xenium_colors)]

xenium_list <- map(rafalib::splitit(Xenium_hBrain$cell_type_class), ~Xenium_hBrain$Genes[.x])

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Xenium_hBrain.pdf")), height = 15, width = 12)
walk(xenium_list,
    ~sce |>
        scDotPlot(features = .x,
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "Xenium",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                    "Xenium" = xenium_colors),
                  clusterRows = FALSE,
                  groupLegends = FALSE
        ) |>
        print())
dev.off()

#### Plot on Visium Data ####

rowData(spe)$Xenium <- NULL
rowData(spe)$Xenium <- Xenium_hBrain$anno[match(rownames(spe), Xenium_hBrain$Genes)] 
table(rowData(spe)$Xenium)

pdf(here(plot_dir, sprintf("SpD_dotplot_Xenium_hBrain.pdf")), height = 14, width = 8)
walk(xenium_list,
     ~spe |>
         scDotPlot(features = .x,
                   group = "SpD",
                   groupAnno = "SpD",
                   featureAnno = "Xenium",
                   scale = TRUE,
                   annoColors = list("SpD" = SpD_colors,
                                     "Xenium" = xenium_colors),
                   clusterRows = FALSE,
                   groupLegends = FALSE
         ) |>
         print())
dev.off()

#### Cell Type Enrichment ####

enrichment_stats_top <- read_csv(here("processed-data", "04_snRNA-seq", "31_sn_subcluster_heatmap", "sce_subcluster_enrichment_top5.csv")) |>
    dplyr::rename(Genes = ensembl)

non_unique <- enrichment_stats_top |> count(ensembl) |> filter(n != 1) |> pull(ensembl)

enrichment_stats_top |> inner_join(Xenium_hBrain)

enrichment_stats_top |> filter(Genes %in% Xenium_hBrain$Genes)

enrichment_stats_select <- enrichment_stats_top |> 
    filter(cell_type_anno %in% c("Astro.1", "Astro.2", "Astro.3",
                                 "Oligo.1", "Oligo.2", "Oligo.3", "Oligo.5",
                                 "OPC.5")
           ) |> 
    dplyr::rename(gene_name = Genes)|>
    inner_join(rowData(sce) |> as.data.frame()) |> 
    filter(gene_name != "ENSG00000272076") ## top gene for several astros

enrichment_stats_select |> dplyr::count(cell_type_anno, gene_type)

enrichment_stats_select |> filter(gene_type == "protein_coding")

rowData(sce)$enrichment <- NULL
rowData(sce)$enrichment <- enrichment_stats_select$cell_type_anno[match(rownames(sce), enrichment_stats_select$gene_name)] 
table(rowData(sce)$enrichment)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_enrichment_select.pdf")), height = 11, width = 8)
sce |>
    scDotPlot(features = enrichment_stats_select$gene_name,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "enrichment",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "enrichment" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 
dev.off()

#### Single Cell MeanRatio ####

load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio","marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats_MeanRatio |>
    filter(MeanRatio.rank <= 10, MeanRatio > 1) |>
    arrange(cellType.target)

## 11 top5 genes overlap
mean_ratio_xenium_overlap <- marker_stats_top |> 
    inner_join(Xenium_hBrain |> 
                   dplyr::rename(gene = Genes)) |>  
    select(gene, cellType.target, MeanRatio, MeanRatio.rank, Annotation, cell_type_anno) |>
    arrange(cellType.target) 

mean_ratio_xenium_overlap |> print(n = 21)

marker_stats_select <- marker_stats_top |> 
    filter(MeanRatio.rank <= 5) |>
    filter(cellType.target %in% c("Astro.1", "Astro.2", "Astro.3",
                                 "Oligo.1", "Oligo.2", "Oligo.3", "Oligo.5",
                                 "OPC.5",
                                 "Excit.L5.2", "Excit.L2_5.1")
    )

marker_stats_select |> dplyr::count(cellType.target)

marker_stats_select 

rowData(sce)$MeanRatio <- NULL
rowData(sce)$MeanRatio <- as.character(marker_stats_select$cellType.target)[match(rownames(sce), marker_stats_select$gene)] 
table(rowData(sce)$MeanRatio)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_MeanRatio_select.pdf")), height = 11, width = 8)
sce |>
    scDotPlot(features = marker_stats_select$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "MeanRatio",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "MeanRatio" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 
dev.off()

#### Cell type specific modeling ####

broad_cell_types <- c("Oligo", "Astro")
names(broad_cell_types) <- broad_cell_types

ct_specific_mod <- map(broad_cell_types, ~readRDS(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", .x, sprintf("modeling_results_subtype-%s.rds", .x))) |>
                           pluck("enrichment") |>
                           pivot_longer(!c("ensembl", "gene"), names_to = "stat") |>
                           mutate(stat = gsub("p_value", "p.value", stat),
                                  stat = gsub("t_stat", "t.stat", stat),) |>
                           separate(stat, into = c("stat", "cellType.target"), sep = "_") |>
                           pivot_wider(names_from = "stat", values_from = "value", names_prefix = "enrich_") |>
                           dplyr::rename(gene_ensembl = ensembl)
                       )

head(ct_specific_mod$Oligo) 

ct_specific_marker_stats <- map2(broad_cell_types, ct_specific_mod, ~get(load(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", .x, sprintf("marker_stats_MeanRatio_%s.Rdata", .x)))) |>
                                 left_join(.y))


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


#### Genes of Interest ####

NE_receptor_genes <- c("ADRA1A","ADRA1B","ADRA1D","ADRA2A","ADRA2B","ADRA2C","ADRB1","ADRB2","ADRB3")
all(NE_receptor_genes %in% rownames(sce))

genes_other <- c("RELN")

genes_of_interest <- c(NE_receptor_genes, genes_other)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_genes_of_interest.pdf")), height = 7, width = 8)
sce |>
    scDotPlot(features = genes_of_interest,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 
dev.off()

pdf(here(plot_dir, sprintf("SpD_dotplot_genes_of_interest.pdf")), height = 7, width = 8)
spe |>
    scDotPlot(features = genes_of_interest,
              group = "SpD",
              groupAnno = "SpD",
              scale = TRUE,
              annoColors = list("cell_type_anno" = SpD_colors),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 
dev.off()

#### Mediation Genes ####

# list.files(here("processed-data","22_Mediation","out-erc_astro"))

mediation_summary <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediation_vs_baseline_summary.tsv.gz"))

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) 

## 20 outcomes in base probes
mediator_outcome |> filter(outcome %in% Xenium_hBrain$Genes)

mediator_outcome |> dplyr::count(med_cl, mediator)
mediator_outcome |> dplyr::count(med_cl, mediator)

head(mediation_summary)

mediator_tab <- mediation_summary |> filter(mediated_n > 0) |> select(mediation_run) |> separate(mediation_run, into = c("cell_type", "ensemblID", "gene_name"), sep = "\\|")

mediator_tab$gene_name %in% Xenium_hBrain$genes

pdf(here(plot_dir, "sn_dotplot_mediator_genes.pdf"))
sce |>
    scDotPlot(features = mediator_tab$gene_name,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = FALSE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 

sce |>
    scDotPlot(features = mediator_tab$gene_name,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    ) 
dev.off()
