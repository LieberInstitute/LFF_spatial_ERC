# March 2024 - Louise Huuki-Myers
# Compile, reformat, and compare marker genes from literature 

library("tidyverse")
library("here")

## Prep directories
plot_dir <- here("plots", "05_spe_correct_cluster", "00_lit_marker_genes_layer")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "05_spe_correct_cluster", "00_lit_marker_genes_layer")
if(!dir.exists(data_dir)) dir.create(data_dir)


#### summarize marker genes from lit recorded in reading ####
layer_markers <- read.csv(here("processed-data", "05_spe_correct_cluster","16_plot_layer_markers", "genes_to_plot.csv")) |>
    select(gene_name = gene, Layer = Layer_Simple, Source)

rna_markers <- read.csv(here("processed-data", "05_spe_correct_cluster","16_plot_layer_markers", "rnascope_layer_markers.csv"))
all(rna_markers$gene %in% layer_markers$gene_name)

lit_markers <- read.csv(here("processed-data", "05_spe_correct_cluster","00_lit_marker_genes_layer", "erc_layer_lit_marker_genes.csv")) |>
    filter(!(gene_name %in% layer_markers$gene_name & grepl("Maynard", source))) |>
    select(gene_name, Layer = layer, Source = source) |>
    mutate(Layer = gsub("L", "Layer", Layer))

lit_marker_summary <- rbind(layer_markers, lit_markers) |>
    group_by(gene_name, Layer) |>
    summarize(n_studies = n(),
              studies = paste0(Source, collapse = ",")) |>
    arrange(Layer)

lit_marker_summary |> ungroup() |> count(Layer)

lit_marker_summary |> write_csv(here(data_dir, "lit_layer_marker_summary.csv"))

