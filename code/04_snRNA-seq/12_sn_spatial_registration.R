## January 2025, Louise Huuki-Myers
## Compare ERC spatial Domains to DLPFC layers

library("spatialLIBD")
library("purrr")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "10_spatial_registration_DLPFC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### get reference layer enrichment statistics ####
layer_modeling_results <- map(c(HumanPilot = "modeling_results", spatialDLPFC = "spatialDLPFC_Visium_modeling_results"), fetch_data)

#### load data ####

## Add spatialDLPFC spatial domain annotations to modeling
dlpfc_anno <- read.csv("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv") |>
    dplyr::filter(bayesSpace == "k09") |>
    mutate(layer_combo2 = gsub(" ", "~", layer_combo2))

dlpfc_colnames <- colnames(layer_modeling_results$spatialDLPFC$enrichment)
pwalk(dlpfc_anno, function(...) dlpfc_colnames <<- gsub(..5, ..3, dlpfc_colnames))
colnames(layer_modeling_results$spatialDLPFC$enrichment) <- dlpfc_colnames
head(layer_modeling_results$spatialDLPFC$enrichment)

#### correlate layer stats ####

cor_layer <- layer_stat_cor(stats = registration_t_stats,
                            modeling_results = layer_mod,
                            model_type = "enrichment",
                            top_n = 100)

cor_layer <- cor_layer[,order(colnames(cor_layer))]

#### Annotate ####

annotate_registered_clusters(
    cor_stats_layer = .x,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
))


#### create registration heatmaps ####


# slurmjobs::job_single('10_spatial_registration_DLPFC', create_shell = TRUE, memory = '25G', command = "Rcript 10_spatial_registration_DLPFC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




