## Louise Huuki-Myers, May 2025
## Compare gilal re-clustering and glial sub-clustering results

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("spatialLIBD")
library("tidyverse")
library("jaffelab")
library("DeconvoBuddies")
library("bluster")
library("ComplexHeatmap")

data_dir <- here("processed-data", "04_snRNA-seq", "27_sn_gila_check")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "27_sn_gila_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

sce <- sce[,sce$cell_type_class == "glia"]
message(Sys.time(), sprintf(" - Subset to %s, ncells = %i", "glia", ncol(sce)))

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Add cluster info ####
message(Sys.time(), " - Add cluster info")

## Glia clustering
(cluster_fn <- here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", "glia", "walktrap_snn_k10_subclusters_glia.Rdata"))
load(cluster_fn, verbose = TRUE)

length(clusters) == ncol(sce)

glia_cluster_info <- read_csv(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", "glia", "ERCsn_subtype_clustering_info-glia_k10_anno.csv"))
glia_cluster_info |> count(sctype)

glia_anno_table <- glia_cluster_info |> select(k10, anno) |> column_to_rownames("k10")
sce$glia_subtype <- glia_anno_table[paste0("k", clusters),]

table(sce$glia_subtype)

## subtype clustering 
glia_cell_types <- c("Astro", "Endo", "Micro", "Oligo", "OPC")
names(glia_cell_types) <- glia_cell_types

sce$glia_subcluster <- NA

subcluster_info <- map_dfr(glia_cell_types, function(ct){
    ## load clusters
    cluster_fn <- here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", ct, sprintf("walktrap_snn_k10_subclusters_%s.Rdata", ct))
    load(cluster_fn)
    
    ## load cluster info
    subcluster_info <- read_csv(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", ct , sprintf("ERCsn_subtype_clustering_info-%s_k10.csv", ct)),
                                show_col_types = FALSE) |>
        mutate(cell_type_broad = ct)
    
    anno_table <- subcluster_info |> select(k10, cell_type_k10) |> column_to_rownames("k10")
    
    sce[,sce$cell_type_broad == ct]$glia_subcluster <<- anno_table[paste0("k", clusters),]

    return(subcluster_info)
})

table(sce$glia_subcluster)

table(sce$glia_subcluster, sce$glia_subtype)

#### DE compat check ####

de_compat <- colData(sce) |>
    as.data.frame() |>
    count(glia_subcluster, sample_id, APOE_carrier) |>
    filter(n >= 10) |> ## 10 or more nuc
    count(glia_subcluster, APOE_carrier) |>
    filter(n >=2) |> ## 2 or more donors
    count(glia_subcluster) |>
    filter(n >=2) |> ## in both APOE groups
    pull(glia_subcluster)

subcluster_info <- subcluster_info |> mutate(DGE_compat = cell_type_k10 %in% de_compat)

subcluster_info |> filter(!DGE_compat)

## calc concordance
message("Concordance rand score: ", bluster::pairwiseRand(sce$glia_subcluster, sce$glia_subtype, mode = "index"))
# 0.48

#### define colors for subtypes ####

subclusters <- unique(sce$glia_subcluster)

subclusters <- sort(subclusters[!grepl("QC", subclusters)])

subclusters_colors <- create_cell_colors(cell_types = subclusters, split = "\\.")

#### plot marker genes for subcluster ####
message(Sys.time(), " - Plot lit marker genes")

lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    rename(cell_type_target = cell_type) |>
    filter(gene_name %in% rowData(sce)$gene_name)

rownames(sce) <- rowData(sce)$gene_name

lit_markers_list <- map(splitit(lit_markers$cell_type_target), ~lit_markers$gene_name[.x])

plot_marker_express_List(sce[,!grepl("QC", sce$glia_subcluster)], 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_sn_lit_markers_subtype_glia_subcluster.pdf"),
                         cellType_col = "glia_subcluster",
                         gene_name_col = "gene_name",
                         color_pal = subclusters_colors$fine
                         )

#### Jaccard matrix heatmap ####
message(Sys.time(), " - Jaccard matrix heatmap")

jacc.mat <- linkClustersMatrix(sce$glia_subcluster, sce$glia_subtype)

## subcluster annotation
subcluster_anno <- subcluster_info |>
    select(cell_type_k10, metric_fail, cell_type_broad, n) |>
    column_to_rownames("cell_type_k10")

identical(rownames(jacc.mat), rownames(subcluster_anno))

subcluster_ha_row <- rowAnnotation(df = subcluster_anno,
                                   col = list(cell_type_broad = cell_type_colors$broad,
                                              metric_fail = c(detected = "firebrick1",
                                                              `detected,sum` = "firebrick", 
                                                              scDblFinder.score = "sienna")))

## subtype annotation
subtype_anno <- glia_cluster_info |>
    mutate(cell_type_broad = ss(anno, "\\.")) |>
    select(anno, metric_fail, cell_type_broad, n) |>
    column_to_rownames("anno")

identical(colnames(jacc.mat), rownames(subtype_anno))

subtype_ha_col <- HeatmapAnnotation(df = subtype_anno,
                                    col = list(cell_type_broad = c(cell_type_colors$broad, VLMC = "maroon"),
                                               metric_fail = c(detected = "firebrick1",
                                                               `detected,sum` = "firebrick",
                                                               scDblFinder.score = "sienna")))

## plot heatmap
pdf(here(plot_dir, "glia_subtype_v_subcluster_jaccard_heatmap.pdf"), width = 10, height = 10)
Heatmap(jacc.mat,                  name = "Correspondence",
        right_annotation = subcluster_ha_row,
        bottom_annotation = subtype_ha_col,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = c("black", viridisLite::plasma(100)),
        na_col = "black")
dev.off()


#### run cell type specific modeling ####
modeling_fn <- here(data_dir, "modeling_results_subtype_subcluster.rds")
pseudobulk_fn = here(data_dir, "sce_pseudobulk_subtype_subcluster.rds")

if(file.exists(modeling_fn) & file.exists(pseudobulk_fn)){
    
    message(Sys.time(), " - Load modeling data")
    modeling_results <- readRDS(modeling_fn)
    
} else {
    
    message(Sys.time(), " - Running Spatial Registration on: subcluster")

    ## Filter to passALL_metricQC nuclei
    sce <- sce[,sce$passALL_metricQC]
    
    ## make APOE syntatic
    sce$APOE <- gsub("/", "", sce$APOE)
    
    modeling_results <-registration_wrapper(
        sce = sce[,!grepl("QC", sce$glia_subcluster)],
        var_registration = "glia_subcluster",
        var_sample_id = "sample_id",
        covars = c("APOE", "Sex", "Age", "Anc_Afr"),
        gene_ensembl = "gene_id",
        gene_name = "gene_name",
        min_ncells = 10,
        pseudobulk_rds_file = pseudobulk_fn
    )
    
    message(Sys.time(), " - Saving Data")
    saveRDS(modeling_results, file = modeling_fn)
    
}

sce_pseudo <- readRDS(pseudobulk_fn)

top_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_genes, here(data_dir, sprintf("subcluster_enrichment_top10.csv")))

#### spatial registration ####

layer_modeling_results <- list()
## sestan sn enrichment
layer_modeling_results$sestan_EC <- readRDS(here("processed-data", "04_snRNA-seq", "24_external_data_check", "sestan_EC_modeling.rds"))
## PsychENCODE enrichment
layer_modeling_results$snDLPFC_PEC <- list()
layer_modeling_results$snDLPFC_PEC$enrichment <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/registration_stats_LIBD.rds")

## Spatial ERC
layer_modeling_results$spatialERC <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
colnames(layer_modeling_results$spatialERC$enrichment) <- gsub("_Sp", "~Sp", colnames(layer_modeling_results$spatialERC$enrichment))
colnames(layer_modeling_results$spatialERC$enrichment) 

#### correlate layer stats ####
cor_layer <- map(layer_modeling_results, function(layer_mod){
    cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                                modeling_results = layer_mod,
                                model_type = "enrichment",
                                top_n = 100)
    return(cor_layer)
})

#### Registration Annotation ####
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

## save data
spatial_registration <- list(cor_layer = cor_layer, anno = anno)
save(spatial_registration, file = here(data_dir, sprintf("ERCsn_subcluster_registration.Rdata")))

# load(here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("ERCsn_subtype_registration-%s.Rdata", cell_type)))

walk(names(cor_layer), function(ref){
    message(ref)
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_subcluster.pdf", ref)), 
        width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        # reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = subclusters_colors$fine
    ))
    dev.off()
})

subcluster_info2 <- subcluster_info |> left_join(anno_summary |> rename(cell_type_k10 = cluster, sestan_EC_c = sestan_EC, snDLPFC_PEC_c = snDLPFC_PEC))

write_csv(subcluster_info2, file = here(data_dir, "ERCsn_glia_subcluster_info.csv"))

subcluster_info2 |> 
    # group_by(cell_type_broad, passALL_metricQC) |> 
    # group_by(passALL_metricQC) |>
    filter(passALL_metricQC) |>
    group_by(cell_type_broad) |>
    summarise(n_cluster = n(),
              n_nuclei = sum(n))

glia_cluster_info |> 
    mutate(cell_type_broad = ss(anno,"\\.")) |>
    # group_by(passALL_metricQC) |> 
    # group_by(cell_type_broad, passALL_metricQC) |> 
    filter(passALL_metricQC) |> 
    group_by(cell_type_broad) |>
    summarise(n_cluster = n(),
              n_nuclei = sum(n))

# slurmjobs::job_single(name = "27_sn_glia_check", memory = "50G", command = "Rscript 27_sn_glia_check.R", create_shell = TRUE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



