## Louise Huuki-Myers, July 2025
## Make dotplots of select gene sets

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("scDotPlot")
library("tidyverse")

plot_dir <- here("plots", "04_snRNA-seq", "10.5_cluster_dotplot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- loadHDF5SummarizedExperiment(dir = here("processed-data", "sce_objects", "sce_reprocess"))
sce

rownames(sce) <- rowData(sce)$Symbol

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors
cell_type_colors$broad <- c(cell_type_colors$broad, Neuron = "lightblue")

#### lit marker dot plot ####

lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    mutate(cell_type_broad = factor(cell_type_broad, levels = names(cell_type_colors$broad))) |>
    filter( gene_name %in% rowData(sce)$Symbol,
            cell_type_broad %in% names(cell_type_colors$broad),
           !is.na(cell_type_broad),
           grepl("Grubman|Davila-Velderrain|Li", studies)
           ) |>
    arrange(cell_type_broad) 

lit_markers |> print(n = 50)

all(lit_markers$gene_name %in% rownames(sce))

lit_markers |> ungroup() |> count(cell_type_broad)
lit_markers |> count(gene_name) |> filter(n > 1)
lit_markers |> filter(gene_name == "CLDN5")

cell_type_colors$broad <- c(cell_type_colors$broad, Neuron = "lightblue")
unique(lit_markers$cell_type_broad)[!unique(lit_markers$cell_type_broad) %in% names(cell_type_colors$broad)]

rowData(sce)$LitMarker <- NULL
rowData(sce)$LitMarker <- lit_markers$cell_type_broad[match(rownames(sce), lit_markers$gene_name)] 
table(rowData(sce)$LitMarker)


pdf(here(plot_dir, sprintf("sn_ct_broad_k10_dotplot_lit_genes.pdf")), height = 9, width = 6)
sce |>
    scDotPlot(features = unique(lit_markers$gene_name),
              group = "ct_broad_k10",
              groupAnno = "ct_broad_k10",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("ct_broad_k10" = cell_type_colors$broad,
                                "LitMarker" = cell_type_colors$broad),
              clusterRows = FALSE,
              groupLegends = FALSE
    )
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

