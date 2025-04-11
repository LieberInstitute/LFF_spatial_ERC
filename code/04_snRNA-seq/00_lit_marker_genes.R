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
    bind_rows(pec_marker_sub_tab) |>
    mutate(source = "PsychENCODE AHBA PFC",
           gene = gsub("\\+", "", gene))


pec_marker_tab |> dplyr::count(cell_type_broad)
# cell_type_broad     n
# <chr>           <int>
# 1 Astro             221
# 2 Endo              287
# 3 Excit             232
# 4 Inhib             241
# 5 Micro             285
# 6 OPC               122
# 7 Oligo             255
# 8 Tcell               2
# 9 VLMC              218

pec_marker_tab |> dplyr::count(cell_type_broad, cell_type_fine)

## format for sctype Database

pec_db_neuron <- pec_marker_tab |> 
    filter(is.na(cell_type_fine) & cell_type_broad %in% c("Excit", "Inhib")) |>
    group_by(cell_type_broad) |>
    summarise(geneSymbolBroad = paste(gene, collapse = ",")) 

pec_db <- pec_marker_tab |> 
    filter(!is.na(cell_type_fine) | !cell_type_broad %in% c("Excit", "Inhib")) |>
    mutate(cellName = ifelse(is.na(cell_type_fine), cell_type_broad, paste0(cell_type_broad,"." ,cell_type_fine))) |> 
    group_by(cellName, cell_type_broad, cell_type_fine) |>
    summarise(geneSymbolmore1 = paste(gene, collapse = ",")) |>
    left_join(pec_db_neuron) |>
    mutate(tissueType = "Brain",
           geneSymbolmore2 = 0,
           geneSymbolmore1 = ifelse(is.na(geneSymbolBroad), 
                                geneSymbolmore1, 
                                paste0(geneSymbolBroad, ",", geneSymbolmore1)) ## Add neuron broad markers to subtypes
               ) |>
    select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

pec_db |> print(n=35)

writexl::write_xlsx(pec_db, path = here(data_dir, "scTypeDB_PEC_PFC_custom.xlsx"))

## test
gs_list_pec <- gene_sets_prepare(here(data_dir, "scTypeDB_PEC_PFC_custom.xlsx"), "Brain")

#### sctype genes ####
library("HGNChelper")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
gs_list <- gene_sets_prepare(db_, "Brain")

sctype_marker_tab <- tibble(source_cell_type = rep(names(gs_list$gs_positive), map_int(gs_list$gs_positive, length)),
                            gene = unlist(gs_list$gs_positive),
                            ) |> left_join(tibble(source_cell_type = c(
                                'Astrocytes',
                                'Endothelial cells',
                                'Microglial cells',
                                'Oligodendrocytes',
                                'Oligodendrocyte precursor cells',
                                'Glutamatergic neurons',
                                'GABAergic neurons',
                                'Mature neurons'), 
                                cell_type_broad  = c(
                                    'Astro',
                                    'Endo',
                                    "Micro",
                                    'Oligo',
                                    'OPC',
                                    'Excit',
                                    'Inhib',
                                    'Neu')
                            )) |> mutate(source = "sctype Brain")

sctype_marker_tab |> dplyr::count(cell_type_broad)



#### summarize marker genes from lit recorded in reading ####
lit_markers <- read.csv(here("processed-data", "04_snRNA-seq","00_lit_marker_genes", "lit_marker_genes.csv"))

cell_type_key <- tibble(cell_type = c("Glutamate receptor", "Fibroblast", "Inhib-subtype", "Pericytes", "choroid plexis"),
       cell_type_broad = c("Excit", "Endo", "Inhib", "Endo", "Endo"))

lit_marker_summary <- lit_markers |>
    filter(note != "exclude") |>
    group_by(gene_name, cell_type) |> 
    summarize(n_studies = n(),
              studies = paste0(source, collapse = ",")) |>
    left_join(cell_type_key) |>
    mutate(cell_type_broad = ifelse(is.na(cell_type_broad), cell_type, cell_type_broad)) |>
    arrange(cell_type_broad, cell_type)

lit_marker_summary |> ungroup() |> dplyr::count(cell_type)

lit_marker_summary |> dplyr::count(gene_name) |> arrange(-n)
lit_marker_summary |> filter(gene_name == 'SYT1')

lit_marker_summary |> write_csv(here(data_dir, "lit_marker_summary.csv"))

#### compile all marker genes ####


