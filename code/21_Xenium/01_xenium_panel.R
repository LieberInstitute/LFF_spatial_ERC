## Louise Huuki-Myers, Jan 2026
## Explore existing pannel, build custom list


#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("scDotPlot")

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

#### Explore Xenium data ####

Xenium_hBrain <- read_csv(here("processed-data", "21_Xenium", "Xenium_hBrain_v1_metadata.csv")) |>
    mutate(in_sce = Genes %in% rownames(sce))

Xenium_hBrain |> filter(!in_sce)
# Genes  Ensembl_ID      Num_Probesets Codewords Annotation in_sce
# <chr>  <chr>                   <dbl>     <dbl> <chr>      <lgl> 
# 1 BTBD11 ENSG00000151136             8         1 Pvalb      FALSE

Xenium_hBrain <- Xenium_hBrain |> filter(in_sce)
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

## dotplot of panel on cell_type_anno

rowData(sce)$Xenium <- NULL
rowData(sce)$Xenium <- Xenium_hBrain$Annotation[match(rownames(sce), Xenium_hBrain$Genes)] 
table(rowData(sce)$Xenium)

pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_Xenium_hBrain.pdf")), height = 15, width = 12)
sce |>
    scDotPlot(features = Xenium_hBrain$Genes,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Xenium",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = TRUE,
              groupLegends = FALSE
    )
dev.off()


#### cell type modeling ####
## Global Enrichment 
enrichment_stats_top <- read_csv(here("processed-data", "04_snRNA-seq", "31_sn_subcluster_heatmap", "sce_subcluster_enrichment_top5.csv")) |>
    dplyr::rename(Genes = ensembl)

non_unique <- enrichment_stats_top |> count(ensembl) |> filter(n != 1) |> pull(ensembl)

enrichment_stats_top |> inner_join(Xenium_hBrain)

enrichment_stats_top |> filter(Genes %in% Xenium_hBrain$Genes)


#### Oligo modeling ####

load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio","marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats_MeanRatio |>
    filter(MeanRatio.rank <= 5, MeanRatio > 1) |>
    arrange(cellType.target)

marker_stats_top |> filter(gene %in% Xenium_hBrain$Genes)
