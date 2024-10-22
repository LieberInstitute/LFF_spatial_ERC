## October 2024, Louise Huuki-Myers
## plot layer markers

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "05_spe_correct_cluster", "16_plot_layer_markers")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "16_plot_layer_markers")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load SPE data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))
# spe <- readRDS(here("code", "06_iSEE_app", "sce_ERC_iSEE.rds"))

logcounts(spe)[1:5,1:5]
rownames(spe) <- rowData(spe)$gene_name

## load layer markers
rna_scope_markers <- read_csv(here(data_dir, "rnascope_layer_markers.csv")) |>
    arrange(Source, Layer)

# gene    Layer  Source             
# <chr>   <chr>  <chr>              
# 1 RELM    Layer1 ERC RNAScope       
# 2 FREM3   Layer3 ERC RNAScope       
# 3 TRABD2A Layer5 ERC RNAScope       
# 4 MBP     WM     ERC RNAScope 

dlpfc_layer_markers <- read_csv(here("processed-data", "00_project_prep", "06_marker_genes", "dlpfc_layers_top100.csv")) |>
    rename(Source = dataset) |>
    mutate(layer = gsub("L([0-9])", "Layer\\1", layer)) |>
    arrange(Source, layer) |>
    mutate(layer_combo = ifelse(is.na(SpD), layer, paste0(layer, " ~ ", SpD)))

dlpfc_layer_markers |>
    filter(top == 1, layer == "Layer1") 
    
dlpfc_layer_markers |> filter(!is.na(SpD))

dlpf_top_markers <- dlpfc_layer_markers |>
    filter(top == 1) |>
    select(gene, Layer = layer_combo, Source)

genes_to_plot <- bind_rows(dlpf_top_markers, rna_scope_markers) |> 
    group_by(gene) |>
    summarize(n = n(), 
              Source = paste(sort(Source), collapse = ", "), 
              Layer = paste(sort(Layer), collapse = ", ")) |>
    mutate(Layer_Simple = gsub(",| .*?$","",Layer)) |>
    arrange(Layer)
    
pmap(genes_to_plot, function(gene, Layer_Simple){
    
    filename = sprintf("spe_erc_%s-%s.pdf", Layer_Simple, gene)
    return(filename)
    
})

#### Plot Layer Markers ####

sample_order <- sort(unique(spe$BrNum))
walk(c("MBP"),
     # walk(c("MBP", "AQP4", "HPCAL1", "KRT17", "MOBP", "PCP4", "SNAP25"),
     ~vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, sprintf("spe_erc-%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order)
     
)
