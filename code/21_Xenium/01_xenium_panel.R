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

#### Plot base pannel in snRNA-seq genes ####

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

#### Plot base pannel on Visium Data ####

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







