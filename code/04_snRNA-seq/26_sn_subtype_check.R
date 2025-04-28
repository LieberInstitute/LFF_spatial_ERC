## Louise Huuki-Myers, April 2025
## Run subtype clustering on sn data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("getopt")
library("spatialLIBD")
library("tidyverse")
# library("scater")
# library("scran")
# library("scry")

# Import command-line parameters
scec <- matrix(
    c("cell_type", "c", "1", "character", "Name of cell type to sub-cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

cell_type <- opt$cell_type

# cell_type = "Endo"

data_dir <- here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", cell_type)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "26_sn_subtype_check", cell_type)
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

stopifnot(cell_type %in% levels(sce$cell_type_broad))

## subset by Broad cell type
sce <- sce[,sce$cell_type_broad == cell_type]
message(Sys.time(), sprintf(" - Subset to %s, ncells = %i", cell_type, ncol(sce)))

## read in cluster data
load(here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", cell_type, sprintf("walktrap_snn_k20_subclusters_%s.Rdata", cell_type)),
     verbose = TRUE) 
# clusters

sce$cell_type_fine <- sprintf("%s.%i", cell_type, clusters)
table(sce$cell_type_fine)

#### run cell type specific modeling ####
## make APOE syntatic
sce$APOE <- gsub("/", "", sce$APOE)

cluster_var <- "cell_type_fine"
message(Sys.time(), " - Running Spatial Registration on: ", cluster_var, "-", cell_type)
stopifnot(cluster_var %in% colnames(colData(sce)))

table(sce[[cluster_var]], sce$sample_id)

modeling_results <-registration_wrapper(
    sce = sce,
    var_registration = cluster_var,
    var_sample_id = "sample_id",
    covars = c("APOE", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", cell_type))
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("modeling_results_subtype-%s.rds", cell_type)))

sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", cell_type)))

top_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s.csv", cell_type)))

#### spatial registration ####

layer_modeling_results <- list()
## sestan sn enrichment
layer_modeling_results$sestan_EC <- readRDS(here("processed-data", "04_snRNA-seq", "24_external_data_check", "sestan_EC_modeling.rds"))

#### correlate layer stats ####
cor_layer <- map(layer_modeling_results, function(layer_mod){
    cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                                modeling_results = layer_mod,
                                model_type = "enrichment",
                                top_n = 100)
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


anno_summary2 <- colData(sce) |> 
    as.data.frame() |> 
    group_by(cluster = cell_type_fine)  |> 
    summarise(n_cells = n(), n_donors = length(unique(sample_id))) |>
    left_join(anno_summary)

write_csv(anno_summary2, file = here(data_dir, sprintf("ERCsn_subtype_registration_anno_summary-%s.csv", cell_type)))

## save data
spatial_registration <- list(cor_layer = cor_layer, anno = anno)
save(spatial_registration, file = here(data_dir, sprintf("ERCsn_subtype_registration-%s.Rdata", cell_type)))

# load(here(data_dir, sprintf("ERCsn_subtype_registration-%s.Rdata", cell_type)))

#### create registration heatmaps ####
## load colors
# load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
# load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
# 
# # cell_type_colors
# layer_colors <- list(HumanPilot = spatialLIBD::libd_layer_colors,
#                      spatialDLPFC = NULL, #TODO add spatial domain colors
#                      spatialERC = SpD_colors,
#                      snDLPFC_PEC = NULL,
#                      sestan_EC = NULL)

# map2(cor_layer, layer_colors, ~all(colnames(.x) %in% names(.y)))

walk(names(cor_layer), function(ref){
    message(ref)
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_subtype-%s.pdf", ref, cell_type)), width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        # reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]]
        # query_colors = cell_type_colors$anno
    ))
    dev.off()
})

# slurmjobs::job_loop(loops = list(cell_type = c("Oligo", "Astro", "Micro", "Endo")), create_shell = TRUE, name = "26_sn_subtype_check", create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



