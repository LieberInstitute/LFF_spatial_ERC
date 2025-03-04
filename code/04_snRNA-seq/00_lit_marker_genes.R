# March 2024 - Louise Huuki-Myers
# Compile, reformat, and compare marker genes from literature 

library("tidyverse")
library("HDF5Array")
library("here")
library("rjson")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "00_lit_marker_genes")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "00_lit_marker_genes")
if(!dir.exists(data_dir)) dir.create(data_dir)


#### prep psychENCODE genes ####
pec_json <- fromJSON(file=here("processed-data", "04_snRNA-seq","13_sctype_final", "AHBA_PFC_filtered.json"))

pec_json$cell_types
pec_json[[1]]$subtypes

map_chr(pec_json$cell_types, "name")
# [1] "Excitatory neuron" "Inhibitory neuron" "Astro"             "Oligo"             "OPC"              
# [6] "Endo"              "Micro/PVM"         "T cells"           "VLMC" 

map_lgl(pec_json$cell_types, "subtypes" %in% names())

map_chr(pec_json$cell_types[[1]]$subtypes$cell_types, "name")
# [1] "L2/3 IT"    "L5 IT"      "L6 IT"      "L4 IT"      "L6 CT"      "L6b"        "L5 ET"      "L5/6 NP"   
# [9] "L6 IT Car3"

map(pec_json$cell_types[[1]]$subtypes$cell_types$markers, "genes")

## summarize marker genes from lit recorded in reading
lit_markers <- read.csv(here("processed-data", "04_snRNA-seq", "lit_marker_genes.csv")) |>
    group_by(gene_name, cell_type) |>
    summarize(n_studies = n(),
              studies = paste0(source, collapse = ",")) |>
    mutate(in_data = gene_name %in% rowData(sce)$Symbol)

lit_markers |> write_csv(here("processed-data", "04_snRNA-seq", "lit_marker_summary.csv"))

