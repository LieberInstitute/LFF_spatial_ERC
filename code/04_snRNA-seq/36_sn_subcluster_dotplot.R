## Louise Huuki-Myers, July 2025
## Make dotplots of select gene sets

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("scDotPlot")
library("tidyverse")

data_dir <- here("processed-data", "04_snRNA-seq", "36_sn_subcluster_dotplot")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "36_sn_subcluster_dotplot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors

#### Sex gene dot plot ####

table(as.list(rowRanges(sce)[["seqnames"]]))
sex_check <- list(male = c("SRY", "RPS4Y1", "RPS4Y2", "DDX3Y", "KDM5D", "UTY", "ZFY", "EIF1AY", "USP9Y", "TSPY1"),
                  female = c("XIST", "TSIX", "KDM6A", "EIF2S3X", "RPS4X"))

sex_check <- map(sex_check, ~.x[.x %in% rownames(sce)])

sex_check_tb <- tibble(gene = unlist(sex_check), 
                       sex_check = unlist(map2(names(sex_check), sex_check, ~rep(.x, length(.y)))))

# Autosomal, Hormone-Responsive
# "FOXP3", "ESR1", "AR", "TSHR", "PRLR"

rowData(sce)$sex_check <- sex_check_tb$sex_check[match(rownames(sce), sex_check_tb$gene)] 
table(rowData(sce)$sex_check)


pdf(here(plot_dir, sprintf("sn_sample_dotplot_Sex_genes.pdf")))
sce |>
    scDotPlot(features = unlist(sex_check),
              group = "sample_id",
              groupAnno = "Sex",
              featureAnno = "sex_check",
              # scale = TRUE,
              annoColors = list("Sex" = sex_colors),
              clusterRows = FALSE,
              groupLegends = FALSE
              )
dev.off()


#### risk gene dot plot ####

AD_risk <- read_csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    filter(symbol %in% rowData(sce)$gene_name) 

nrow(AD_risk)


pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_AD_risk.pdf")), height = 10)
sce |>
    scDotPlot(features = unique(AD_risk$symbol),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = TRUE,
              groupLegends = FALSE
    )
dev.off()

#### lit marker dot plot ####

lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    mutate(cell_type_broad = factor(gsub("Endo", "Vasc", cell_type_broad), 
                                    levels = c(levels(sce$cell_type_broad), "Neuron"))) |>
    filter(gene_name %in% rowData(sce)$gene_name,
           !is.na(cell_type_broad)) |>
    # mutate(cell_type_broad = gsub("Endo", "Vasc", cell_type_broad)) |>
    arrange(cell_type_broad) 


lit_markers |> ungroup() |> count(cell_type_broad)
lit_markers |> count(gene_name) |> filter(n > 1)
lit_markers |> filter(gene_name == "CLDN5")

cell_type_colors$broad <- c(cell_type_colors$broad, Neuron = "lightblue")
levels(lit_markers$cell_type_broad)[!unique(lit_markers$cell_type_broad) %in% names(cell_type_colors$broad)]

rowData(sce)$LitMarker <- NULL
rowData(sce)$LitMarker <- lit_markers$cell_type_broad[match(rownames(sce), lit_markers$gene_name)] 
table(rowData(sce)$LitMarker)


pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_lit_genes.pdf")), height = 10, width = 12)
sce |>
    scDotPlot(features = unique(lit_markers$gene_name),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "LitMarker" = cell_type_colors$broad),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()


#### Global Enrichment Dot plot ####
enrichment_stats_top <- read.csv(here("processed-data", "04_snRNA-seq", "31_sn_subcluster_heatmap", "sce_subcluster_enrichment_top5.csv")) |> filter(top == 1)

non_unique <- enrichment_stats_top |> count(ensembl) |> filter(n != 1) |> pull(ensembl)

enrichment_stats_top |> filter(ensembl %in% non_unique) |> count(cell_type_anno)
enrichment_stats_top |> filter(cell_type_anno == "Oligo.03")

enrichment_stats_top_unique <- enrichment_stats_top |> filter(!ensembl %in% non_unique)

rowData(sce)$GlobalEnrich <- NULL
rowData(sce)$GlobalEnrich <- enrichment_stats_top_unique$cell_type_anno[match(rownames(sce), enrichment_stats_top_unique$ensembl)] 
table(rowData(sce)$GlobalEnrich)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Enrich_genes.pdf")), height = 14, width = 10)
sce |>
    scDotPlot(features = enrichment_stats_top_unique$ensembl,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "GlobalEnrich",
              scale = FALSE,
              annoColors = list(cell_type_anno = cell_type_colors$anno,
                                GlobalEnrich = cell_type_colors$anno),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()


#### Mean Ratio Dot plot ####

load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio","marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats_MeanRatio |>
    filter(MeanRatio.rank <= 2, MeanRatio > 1, gene_name %in% rownames(sce_pb)) |>
    arrange(cellType.target)

marker_stats_top  |> print(n = 35)

rowData(sce)$MeanRatio <- NULL
rowData(sce)$MeanRatio <- marker_stats_top$cellType.target[match(rownames(sce), marker_stats_top$gene)] 
table(rowData(sce)$MeanRatio)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_MeanRatio_genes.pdf")), height = 14, width = 10)
sce |>
    scDotPlot(features = marker_stats_top$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "MeanRatio",
              scale = TRUE,
              annoColors = list(cell_type_anno = cell_type_colors$anno,
                                MeanRatio = cell_type_colors$anno),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()


#### Oligo Enrich genes ####
pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Enrich_genes_Oligo.3.pdf")), height = 5, width = 10)
sce |>
    scDotPlot(features = enrichment_stats_top |> filter(cell_type_anno == "Oligo.03") |> pull(ensembl),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "GlobalEnrich",
              scale = TRUE,
              annoColors = list(cell_type_anno = cell_type_colors$anno,
                                GlobalEnrich = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()


#### MeanRatio dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_MeanRatio_genes$cellType.target[match(rownames(sce), enrichment_stats_top$ensembl)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_MeanRatio.pdf", celltype)))
sce |>
    scDotPlot(features = top_MeanRatio_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "Marker" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

#### Oligo OPC marker genes ####

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

## from https://www.biocompare.com/Editorial-Articles/590587-A-Guide-to-Oligodendrocyte-Markers/
oligo_lit_markers <- list(OPC = c("PDGFRA", "CSPG4", "MAG", "CNP", "A2B5"),
                    Oligo = c("PLP1", "ZFP191", "ZFP488", "ZFP536", "SOX17", "NKX6-2", "SMARCA4", "CD82", "TFR", "MAL"),
                    premyelin_Oligo = c("SOX10", "OLIG1", "OLIG2", "NKX2-2", "CD9"),
                    myelinating_Oligo = c("BMP4", "ENPP4", "ASAP", "TMEM10", "MOG"),
                    disease_associated = c("SERPINA3", "C4B", "TNFRSF1A", "IL1B", "IL33", "HMOX1", "TNF", "ERK", "ERK2"), #https://doi.org/10.1038/s41593-025-01873-x
                    AD_risk = c("APOE", "APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1")
)

# MBP, BACE1, and APP capable of producing ABeta  Bright/Ranjani
oligo_lit_markers <- map(oligo_lit_markers, ~.x[.x %in% rownames(sce)])

oligo_lit_markers <- AnnotationDbi::unlist2(oligo_lit_markers)


rowData(sce)$oligo_marker <- NULL
rowData(sce)$oligo_marker <- names(oligo_lit_markers)[match(rownames(sce), oligo_lit_markers)] 
table(rowData(sce)$oligo_marker)

# oligo_lit_markers2 <- list(OPC = c("PDGFRA", "MEGF11"),
#                     Oligo = c("OLIG1", "OLIG2"),
#                     premyelin_Oligo = c("SOX10",  "NKX2-2"),
#                     myelinating_Oligo = c("MBP", "PLP"),
#                     disease_associated = c("SERPINA3", "C4B"), #https://doi.org/10.1038/s41593-025-01873-x
#                     AD_risk = c("APOE", "APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1")
# )


pdf(here(plot_dir, "sn_subtype_OligoOPC5_dotplot_lit_alt_colors.pdf"))
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_anno == "OPC.5"] |>
    scDotPlot(features = oligo_lit_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_marker",
              scale = TRUE,
              # annoColors = list("cell_type_anno" = cell_type_colors$anno),
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## Oligo markers from 
## 10.3389/fcell.2021.714169
oligo_lineage_markers <- list(OPC = c("PDGFRA", "CSPG4"),
                              Oligo_premyelin = c("BMP4", "GPR17"),
                              Oligo_new = c("TMEM2","PROM1"),
                              Oligo = c("BCAS1", "ENPP6", "GJC2"),
                              Oligo_myelin = c("PLP1", "MAG", "MOBP")
)

oligo_lineage_markers <- map(oligo_lineage_markers, ~.x[.x %in% rownames(sce)])

oligo_lineage_markers <- AnnotationDbi::unlist2(oligo_lineage_markers)


rowData(sce)$oligo_lineage_markers <- NULL
rowData(sce)$oligo_lineage_markers <- names(oligo_lineage_markers)[match(rownames(sce), oligo_lineage_markers)] 
table(rowData(sce)$oligo_lineage_markers)


pdf(here(plot_dir, "sn_subtype_OligoOPC5_lineage_dotplot_lit_alt_colors.pdf"))
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_anno == "OPC.5"] |>
    scDotPlot(features = oligo_lineage_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_lineage_markers",
              scale = TRUE,
              # annoColors = list("cell_type_anno" = cell_type_colors$anno),
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

pdf(here(plot_dir, "sn_subtype_OligoOPC5_lineage_dotplot_lit_alt_colors_exp.pdf"))
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_anno == "OPC.5"] |>
    scDotPlot(features = oligo_lineage_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_lineage_markers",
              scale = FALSE,
              # annoColors = list("cell_type_anno" = cell_type_colors$anno),
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

pdf(here(plot_dir, "sn_subtype_OligoOPC_lineage_dotplot_lit_alt_colors.pdf"))
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC"] |>
    scDotPlot(features = oligo_lineage_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_lineage_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

pdf(here(plot_dir, "sn_subtype_OligoOPC_lineage_dotplot_lit_alt_colors_exp.pdf"))
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC"] |>
    scDotPlot(features = oligo_lineage_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_lineage_markers",
              scale = FALSE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

oligo_explore_genes <- c("FGFR3", "PLPP3", "NTSR2", "RYR3") ## Genes of interest from Bernie Jun 26, 2026
all(oligo_explore_genes %in% rownames(sce))

pdf(here(plot_dir, "sn_subtype_OligoOPC_dotplot_explore.pdf"), width = 5)
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC"] |>
    scDotPlot(features = oligo_explore_genes,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "oligo_lineage_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## Oligo Enrichment genes w/ OPC.5 outgroup

oligo_enrichment_genes <- read_csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_enrichment_top10_Oligo.csv")) |>
    select(-`...1`)

rowData(sce)$oligo_enrich_marker <- NULL
rowData(sce)$oligo_enrich_marker <- oligo_enrichment_genes$test[match(rownames(sce), oligo_enrichment_genes$gene)] 
table(rowData(sce)$oligo_enrich_marker)

pdf(here(plot_dir, "sn_subtype_OligoOPC5_dotplot_enrichment_alt_colors.pdf"), width = 5)
# pdf(here(plot_dir, "sn_subtype_OligoOPC5_dotplot_enrichment.pdf"), width = 5)
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC.5"]  |>
    scDotPlot(features = oligo_enrichment_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_enrich_marker",
              scale = TRUE,
              # annoColors = list("cell_type_anno" = cell_type_colors$anno,
              #                   "oligo_enrich_marker" = cell_type_colors$anno),
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "oligo_enrich_marker" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## Oligo markers in ALL cell types
pdf(here(plot_dir, "sn_subtype_dotplot_OligoEnrichment.pdf"), width = 12, height = 9)
sce |>
    scDotPlot(features = oligo_enrichment_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_enrich_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "oligo_enrich_marker" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## top 5 genes

oligo_enrichment_genes |> filter(top <=5) |> pull(gene)

pdf(here(plot_dir, "sn_subtype_OligoOPC5_dotplot_enrichment5_alt_colors.pdf"), width = 5, height =5)

sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_anno == "OPC.5"]  |>
    scDotPlot(features = oligo_enrichment_genes |> filter(top <=5) |> pull(gene),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_enrich_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors,
                                "oligo_enrich_marker" = Oligo_OPC_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

#### Oligo OPC gene list dotplots ####

## blanchard data
source(here("external-data", "Blanchard2022", "get_Blanchard_gene_list.R"))

Blanchard_gene_list <- map(Blanchard_gene_list, ~.x[.x %in% rownames(sce)])
map_int(Blanchard_gene_list, length)


plot_OligoOPC_dotplot <- function(gene_list, name){
    
    pdf(here(plot_dir, sprintf("sn_subtype_OligoOPC5_dotplot_%s_alt_colors.pdf", name)), width = 5, height = (length(gene_list)/3) + 1)
    
    print(sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_anno == "OPC.5"] |>
        scDotPlot(features = gene_list,
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  # featureAnno = "oligo_marker",
                  scale = TRUE,
                  # annoColors = list("cell_type_anno" = cell_type_colors$anno),
                  annoColors = list("cell_type_anno" = Oligo_OPC_colors),
                  clusterRows = TRUE,
                  groupLegends = FALSE))
    
    dev.off()
}

map2(Blanchard_gene_list, names(Blanchard_gene_list), ~plot_OligoOPC_dotplot(.x, .y))


#### Inhib sub-types ####

lit_markers |> filter(cell_type_broad == "Inhib")

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Inhib_logCounts.pdf")), height = 5, width = 5)
sce[, sce$cell_type_broad == "Inhib"]|>
    scDotPlot(features = c(unique(lit_markers |> dplyr::filter(cell_type_broad == "Inhib") |> pull(gene_name)), "PAX6"),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "LitMarker",
              scale = FALSE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
               # clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()

sce[, sce$cell_type_broad == "Inhib"]|>
    scDotPlot(features = c(unique(lit_markers |> dplyr::filter(cell_type_broad == "Inhib") |> pull(gene_name)), "PAX6"),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
               # clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()

#### Micro sub-types ####

lit_markers |> filter(cell_type_broad %in% c("Micro", "Macro"))

# 10.15252/embj.2019101997
# micro_genes <- c("IBA1", "CD11B", "CX3CR1", "CR3", "CD11B", "CX3CR1-EGFP", "CD45", "CCR2", "TREM2", "4D4")

# 10.1083/jcb.201709069
# micro_gene <- c("APOE", "CLU", "SORL1", "ABCA7", "TREM2", "CD33", "MS4A6A", "CR1", "EPHA1", "HLA-DRB1", "IL1RAP")

## 
micro_lit_markers <- c("C1QB", "CSF1R", "CX3CR1", "HLA-DRA", "APOE", "TREM2", "CD33")
micro_lit_markers2 <- micro_lit_markers[micro_lit_markers %in% rownames(sce)]

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Micro.pdf")), height = 5, width = 5)
sce[, sce$cell_type_broad == "Micro" | sce$cell_type_broad == "Macro"]|>
    scDotPlot(features = micro_lit_markers2,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Micro_logCounts.pdf")), height = 5, width = 5)
sce[, sce$cell_type_broad == "Micro" | sce$cell_type_broad == "Macro"]|>
    scDotPlot(features = micro_lit_markers2,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "LitMarker",
              scale = FALSE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()

#### Astro sub-types ####

lit_markers |> filter(cell_type_broad == "Astro")

Astro_lit_markers <- list(disease_associated = c("SERPINA3", "TNFRSF1A", "IL1B", "IL33", "HMOX1"), #https://doi.org/10.1038/s41593-025-01873-x
                          AD_risk = c("APOE", "APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1"),
                          Protoplasmic = c("SLC1A2", "SLC1A3"), #SLC from Dai et al. https://doi.org/10.1186/s40478-023-01624-8
                          Fibrous = c("GFAP", "CD44", "AQP4", "CPAMD8"),
                          Reactive = c("VIM"),
                          Homeostatic = c("NRXN1"),
                          Siletti_cortical = c("WIF1",   # grey matter 
                                               "TNC",    # white matter 
                                               "LMO2") # Interlaminar 
)

Astro_lit_markers <- AnnotationDbi::unlist2(Astro_lit_markers)

rowData(sce)$astro_marker <- NULL
rowData(sce)$astro_marker <- names(Astro_lit_markers)[match(rownames(sce), Astro_lit_markers)] 
table(rowData(sce)$astro_marker)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Astro_logCount.pdf")), height = 5, width = 5)
sce[, sce$cell_type_broad == "Astro"]|>
    scDotPlot(features = Astro_lit_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "astro_marker",
              scale = FALSE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Astro.pdf")), height = 5, width = 5)
sce[, sce$cell_type_broad == "Astro"]|>
    scDotPlot(features = Astro_lit_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "astro_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()

Astro_explore_genes <- read_csv(here(data_dir, "astrocyte_marker_genes_claude.csv")) |>
    select(-ERC_specific_notes) |>
    filter(!gene %in% rownames(sce))

Astro_explore_genes |> filter(gene %in% Astro_lit_markers)

Astro_explore_genes |> filter()

# rowData(sce)$astro_marker <- NULL
# rowData(sce)$astro_marker <- names(Astro_lit_markers)[match(rownames(sce), Astro_lit_markers)] 
# table(rowData(sce)$astro_marker)


test <- sce[, sce$cell_type_broad == "Astro"]|> 
    scDotPlot(features = Astro_explore_genes|> filter(contamination_risk == "none") |> pull(gene),
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "astro_marker",
                  scale = TRUE,
                  annoColors = list("cell_type_anno" = cell_type_colors$anno),
                  clusterRows = TRUE,
                  groupLegends = FALSE
)

p_build <- ggplot_build(test)


pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Astro_explore.pdf")))

sce[, sce$cell_type_broad == "Astro"]|>
    scDotPlot(features = Astro_explore_genes|> filter(contamination_risk == "none") |> pull(gene),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "astro_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce[, sce$cell_type_broad == "Astro"]|>
    scDotPlot(features = Astro_explore_genes |> filter(contamination_risk != "none") |> pull(gene),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "astro_marker",
              scale = FALSE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE
    )

dev.off()

# Allen brain atlas Astros 
# WIF1+ — grey matter astrocytes
# TNC+ — white matter astrocytes
# LMO2+ — interlaminar astrocytes

# Bhattacharjee et al. Layer 1 Astros superficial astrocyte state enriched for Id1 and Id3

astro_genes2 <- c("GFAP", "ITGB4", "ANGPTL4", "SLC1A2", "WIF1", "TNC", "LMO2",
                  "ID1", "ID3")
all(astro_genes2 %in% rownames(sce))
astro_genes2 %in% Astro_explore_genes$gene

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Astro_explore2.pdf")))

sce[, sce$cell_type_broad == "Astro"]|>
    scDotPlot(features = astro_genes2[astro_genes2 %in% rownames(sce)],
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "astro_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()


