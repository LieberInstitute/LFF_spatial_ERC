## Louise Huuki-Myers, Dec 2025
## Plot gene expression of other Oligo clusters

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
# library("scDotPlot")

data_dir <- here("processed-data", "19_other_Oligo", "03_other_Oligo_explore")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "03_other_Oligo_explore")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
scec <- matrix(
    c("dataset", "d", "1", "character", "dataset",
      "cluster", "c", "2", "character", "clustering data",
      "opc", "o", "4", "logical", "include OPCs"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt <- list()
# opt$dataset <- "spatialHPC"
# opt$dataset <- "spatialDLPFC"
# opt$dataset <- "spatialdACC"
# opt$cluster <- "k10"
# opt$opc = TRUE

opt$opc <- as.logical(opt$opc)

if(opt$opc){
    opt$dataset <- paste0(opt$dataset, "_wOPC")
}

message(Sys.time(), " - Data: ",  opt$dataset, ", Cluster: ", opt$cluster, " Include OPC: ", opt$opc)

#### load modeling data ####

modeling_fn <- here("processed-data", "19_other_Oligo", "02_other_Oligo_model", sprintf("modeling_results_Oligo_subtype-%s_%s.rds", opt$dataset, opt$cluster))
pseudobulk_fn = here("processed-data", "19_other_Oligo", "02_other_Oligo_model", sprintf("sce_pseudobulk_Oligo_subtype-%s_%s.rds", opt$dataset, opt$cluster))

message(Sys.time(), " - Load modeling data")
modeling_results <- readRDS(modeling_fn)
    
#### pairwise registration ###

modeling_fn <- list.files(here("processed-data", "19_other_Oligo", "02_other_Oligo_model"), pattern = "modeling_results")
modeling_wOPC_fn <- modeling_fn[grepl("wOPC", modeling_fn)]

names(modeling_wOPC_fn) <- gsub("modeling_results_Oligo_subtype-|.rds", "", modeling_wOPC_fn)

modeling_wOPC <- map(modeling_wOPC_fn, ~readRDS(here("processed-data", "19_other_Oligo", "02_other_Oligo_model", .x)))

cor_layer_pw <- layer_stat_cor(stats = modeling_wOPC$spatialDLPFC_wOPC_k10$enrichment,
                                     modeling_results = modeling_wOPC$spatialdACC_wOPC_k20,
                                     model_type = "enrichment",
                                     top_n = 100)

anno_pw <- annotate_registered_clusters(
    cor_stats_layer = cor_layer_pw,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
)

pdf(here(plot_dir, "other_OligoOPC_pw_layer_stat_cor.pdf"))
layer_stat_cor_plot(cor_layer_pw,
                    # query_colors = other_oligo_colors,
                    # reference_colors = Oligo_OPC_colors2,
                    annotation = anno_pw
)
dev.off()

#### spatialHPC mod ####

# library("ExperimentHub")
# ehub <- ExperimentHub()
# ## Load the HPC dataset
# myfiles <- query(ehub, "humanHippocampus2024")
# hpc_spe <- sce <- myfiles[["EH9605"]]

## spatailHPC modeling
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/08_pseudobulk/PRECAST/visiumHE_DE_stats_domain.rda", verbose = TRUE)
# stats

colnames(stats$enrichment)

cor_layer_spatialHPC <- layer_stat_cor(stats = stats$enrichment,
                               modeling_results = modeling_results,
                               model_type = "enrichment",
                               top_n = 100)

anno_spatialHPC <- annotate_registered_clusters(
    cor_stats_layer = cor_layer_spatialHPC,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
)

pdf(here(plot_dir, "other_Oligo_spatialHPC_layer_stat_cor.pdf"))
layer_stat_cor_plot(cor_layer_spatialHPC,
                    # query_colors = other_oligo_colors,
                    # reference_colors = Oligo_OPC_colors2,
                    annotation = anno_spatialHPC
)
dev.off()

#### Explore proportions ####

cluster_tab_fn <- list(spatialdACC = "walktrap_snn_k20_subclusters_spatialdACC_wOPC.Rdata",
                       spatialDLPFC = "walktrap_snn_k10_subclusters_spatialDLPFC_wOPC.Rdata",
                       spatialHPC = "walktrap_snn_k10_subclusters_spatialHPC.Rdata")

cluster_tab <- map_dfr(cluster_tab_fn, ~get(load(here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster", .x)))) |> 
    mutate(cluster_anno2 = gsub("OPC.[0-9]","OPC", cluster_anno))

## get Oligo annotations
cor_layer_anno_fn <- list(spatialdACC = "other_OligoOPC_cor_layer_anno_spatialdACC_wOPC_k20.Rds",
                          spatialDLPFC = "other_OligoOPC_cor_layer_anno_spatialDLPFC_wOPC_k10.Rds",
                          spatialHPC = "other_Oligo_cor_layer_anno_spatialHPC_k10.Rds")

cor_layer_anno <- map_dfr(cor_layer_anno_fn, ~pluck(readRDS(here("processed-data", "19_other_Oligo", "02_other_Oligo_model", .x)), 'anno')) |>
    dplyr::rename(cluster_anno2 = cluster, anno_confidence = layer_confidence, ERC_Oligo = layer_label)

cluster_sum <- cluster_tab  |>
    dplyr::count(cluster_anno2) |>
    left_join(cor_layer_anno) 

erc_oligo_summary <- read.csv(here("processed-data", "04_snRNA-seq", "33_sn_subcluster_summary", "ERC_sn_subcluster_summary_cell_type.csv")) |>
    filter(cell_type_broad == "Oligo") |>
    select(oligo_anno = cell_type_anno, n = n_nuclei) |>
    mutate(prop = n/sum(n), 
           dataset= "erc", 
           ERC_Oligo = oligo_anno)

cluster_summary <- cluster_sum |>
    filter(!grepl("OPC", cluster_anno2)) |>
    separate(cluster_anno2, sep = "_", into = c("dataset", "oligo_anno"), remove = FALSE) |>
    group_by(dataset) |>
    mutate(prop = n/sum(n)) |>
    bind_rows(erc_oligo_summary)

erc_cluster_summary <- cluster_summary |> 
    ungroup() |>
    group_by(dataset, ERC_Oligo) |>
    summarise(n = sum(n),
              prop = sum(prop),
              clus = paste(c(cluster_anno2), collapse = ", "))

erc_cluster_summary |> print(n = 26)

write_csv(erc_cluster_summary, here(data_dir, "other_Oligo_erc_annotation_summary.csv"))

(erc_cluster_prop_wide <- erc_cluster_summary |> 
    select(-n, -clus) |> 
    pivot_wider(names_from = dataset, values_from = "prop")|>
    arrange(ERC_Oligo))

write_csv(erc_cluster_prop_wide, here(data_dir, "other_Oligo_erc_prop_wide.csv"))


erc_cluster_summary |> select(-prop) |> pivot_wider(names_from = dataset, values_from = "n")

erc_cluster_summary |> print(n = 22)


                   