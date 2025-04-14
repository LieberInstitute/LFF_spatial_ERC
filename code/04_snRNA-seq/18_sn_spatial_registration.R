## January 2025, Louise Huuki-Myers
## Compare ERC cell type populations to DLPFC & ERC spatial domains

library("spatialLIBD")
library("purrr")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "04_snRNA-seq", "18_sn_spatial_registration_anno")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "18_sn_spatial_registration_anno")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load sn ERC modeling ####
erc_sn_modeling_results <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk", "modeling_results-cell_type_anno.rds"))

cell_types <- colnames(erc_sn_modeling_results$enrichment)[grepl("t_stat", colnames(erc_sn_modeling_results$enrichment))]
cell_types <- gsub("t_stat_", "", cell_types)

#### get reference layer enrichment statistics ####
layer_modeling_results <- map(c(HumanPilot = "modeling_results", 
                                spatialDLPFC = "spatialDLPFC_Visium_modeling_results"), 
                              fetch_data)

layer_modeling_results$spatialERC <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
colnames(layer_modeling_results$spatialERC$enrichment) <- gsub("_Sp", "~Sp", colnames(layer_modeling_results$spatialERC$enrichment))
colnames(layer_modeling_results$spatialERC$enrichment) 

## Add spatialDLPFC spatial domain annotations to spatialDLPFC modeling
dlpfc_anno <- read.csv("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv") |>
    dplyr::filter(bayesSpace == "k09") |>
    mutate(layer_combo2 = gsub(" ", "~", layer_combo2))

dlpfc_colnames <- colnames(layer_modeling_results$spatialDLPFC$enrichment)
pwalk(dlpfc_anno, function(...) dlpfc_colnames <<- gsub(..5, ..3, dlpfc_colnames))
colnames(layer_modeling_results$spatialDLPFC$enrichment) <- dlpfc_colnames
# head(layer_modeling_results$spatialDLPFC$enrichment)

## PsychENCODE enrichment
layer_modeling_results$snDLPFC_PEC <- list()
layer_modeling_results$snDLPFC_PEC$enrichment <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/registration_stats_LIBD.rds")

#### correlate layer stats ####
cor_layer <- map(layer_modeling_results, function(layer_mod){
    cor_layer <- layer_stat_cor(stats = erc_sn_modeling_results$enrichment,
                                modeling_results = layer_mod,
                                model_type = "enrichment",
                                top_n = 100)
    
    cor_layer <- cor_layer[cell_types, order(colnames(cor_layer))]
    
    return(cor_layer)
})

#### Annotate ####
anno <- map(cor_layer, ~annotate_registered_clusters(
    cor_stats_layer = .x,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
))

anno_summary <- map2_dfr(anno, names(anno), ~.x |> 
                             select(cluster, layer_label) |>
                             mutate(dataset = .y)) |>
    pivot_wider(names_from = "dataset",
                values_from = "layer_label") |>
    arrange(cluster)

write_csv(anno_summary, file = here(data_dir, "ERCsn_spatial_registration_anno_summary.csv"))

## save data
spatial_registration <- list(cor_layer = cor_layer, anno = anno)
save(spatial_registration, file = here(data_dir, "ERCsn_spatial_registration.Rdata"))

#### create registration heatmaps ####
## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

# cell_type_colors
layer_colors <- list(HumanPilot = spatialLIBD::libd_layer_colors,
                  spatialDLPFC = NULL, #TODO add spatial domain colors
                  spatialERC = SpD_colors,
                  snDLPFC_PEC = NULL)


map2(cor_layer, layer_colors, ~all(colnames(.x) %in% names(.y)))

map(names(cor_layer), function(ref){
    message(ref)
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s.pdf", ref)))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = cell_type_colors$anno
    ))
    dev.off()
})

## dont cluster cols - order by layer
# cell_type_colors$anno <- cell_type_colors$anno[names(cell_type_colors$anno) %in% rownames(cor_layer$spatialERC)]
# cor_layer$spatialERC <- cor_layer$spatialERC[names(cell_type_colors$anno),names(SpD_colors)]
# 
# cor_layer$spatialERC[,names(SpD_colors)]
# cor_layer$spatialERC[names(cell_type_colors$anno),]

names(SpD_colors) %in% colnames(cor_layer$spatialERC)

map(names(cor_layer), function(ref){
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_uncluster.pdf", ref)))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = cell_type_colors$anno,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
    ))
    dev.off()
})


# slurmjobs::job_single('18_sn_spatial_registration', create_shell = TRUE, memory = '25G', command = "Rscript 18_sn_spatial_registration.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




