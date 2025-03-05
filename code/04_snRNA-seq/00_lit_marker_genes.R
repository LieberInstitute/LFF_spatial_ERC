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
pec_json <- fromJSON(file=here("processed-data", "04_snRNA-seq","00_lit_marker_genes", "AHBA_PFC_filtered.json"))

pec_json$cell_types
pec_json$cell_types[[1]]$subtypes

(pec_broad_cell_types <- map_chr(pec_json$cell_types, "name"))
# [1] "Excitatory neuron" "Inhibitory neuron" "Astro"             "Oligo"             "OPC"              
# [6] "Endo"              "Micro/PVM"         "T cells"           "VLMC" 

## broad marker genes
pec_broad_markers <- map(pec_json$cell_types, ~.x$markers[[1]]$genes)

pec_marker_broad_tab <- tibble(source_cell_type = rep(pec_broad_cell_types, map_int(pec_broad_markers, length)),
                         gene = unlist(pec_broad_markers))

## subtype markers
pec_marker_sub_tab <- map2_dfr(c(Excit = 1, Inhib = 2), c("Excit", "Inhib"), function(i, ct){
    subtypes <- map_chr(pec_json$cell_types[[i]]$subtypes$cell_types, "name")
    
    sub_markers <- map(pec_json$cell_types[[i]]$subtypes$cell_types, ~.x$markers[[1]]$genes)
    
    sub_marker_tab <- tibble(source_cell_type = rep(subtypes, map_int(sub_markers, length)),
                             cell_type_broad = ct,
                             cell_type_fine = gsub(" |/", "_", source_cell_type),
                             gene = unlist(sub_markers))
    return(sub_marker_tab)
})

pec_marker_broad_tab |> dplyr::count(source_cell_type)
pec_marker_sub_tab |> dplyr::count(source_cell_type)

pec_marker_tab <- pec_marker_broad_tab |>
    left_join(tibble(source_cell_type = c("Excitatory neuron",
                                          "Inhibitory neuron",
                                          "Micro/PVM",
                                          "T cells"),
                     cell_type_broad = c("Excit",
                                         "Inhib",
                                         "Micro",
                                         "Tcell"))) |>
    mutate(cell_type_broad = ifelse(is.na(cell_type_broad), source_cell_type, cell_type_broad)) |>
    bind_rows(pec_marker_sub_tab)


pec_marker_tab |> dplyr::count(cell_type_broad, cell_type_fine)

## Inhib
map_chr(pec_json$cell_types[[2]]$subtypes$cell_types, "name")
# [1] "Lamp5"      "Lamp5 Lhx6" "Pvalb"      "Chandelier" "Pax6"       "Sncg"       "Sst Chodl"  "Sst"        "Vip"   


#### sctype genes ####

## summarize marker genes from lit recorded in reading
lit_markers <- read.csv(here("processed-data", "04_snRNA-seq", "lit_marker_genes.csv")) |>
    group_by(gene_name, cell_type) |>
    summarize(n_studies = n(),
              studies = paste0(source, collapse = ",")) |>
    mutate(in_data = gene_name %in% rowData(sce)$Symbol)

lit_markers |> write_csv(here("processed-data", "04_snRNA-seq", "lit_marker_summary.csv"))

