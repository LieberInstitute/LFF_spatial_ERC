## October 2024, Louise Huuki-Myers
## Extract layer marker genes from previous studies in easy to access + annotated table

library("spatialLIBD") 
library("tidyverse")
library("here")

#### Set up dirs ####
data_dir <- here("processed-data", "00_project_prep", "06_marker_genes")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


#### Human Pilot manual annotations ####
pilot_layer_top100 <- sig_genes_extract(
    n = 100,
    modeling_results = fetch_data(type = "modeling_results"),
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = fetch_data(type = "sce_layer")
)

dim(pilot_layer_top100)


#### SpatialDLPFC k9 ####
spatialDLPFC_k9_top100 <- sig_genes_extract(
    n = 100,
    modeling_results = fetch_data(type = "spatialDLPFC_Visium_modeling_results"),
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = fetch_data(type = "spatialDLPFC_Visium_pseudobulk")
) 

dlpfc_anno <- read.csv(here("processed-data", "00_project_prep", "spatialDLPFC_Data", "bayesSpace_layer_annotations.csv")) |>
    dplyr::filter(bayesSpace == "k09") |>
    select(SpD = cluster,
           layer = layer_annotation,
           layer_combo)

spatialDLPFC_k9_top100_anno <- spatialDLPFC_k9_top100 |>
    rename(SpD = test) |>
    left_join(dlpfc_anno) |>
    mutate(dataset = "spatialDLPFC",
           marker_anno = layer_combo)

dlpfc_layers_top100 <- pilot_layer_top100 |>
    mutate(dataset = "HumanPilot",
           marker_anno = test) |>
    select(ensembl, gene, dataset, layer=test, top, marker_anno) |>
    bind_rows(spatialDLPFC_k9_top100_anno |>
                  select(ensembl, gene, dataset, SpD, layer, layer_combo, top, marker_anno)) |>
    as_tibble()

write_csv(dlpfc_layers_top100, file = here(data_dir, "dlpfc_layers_top100.csv"))

