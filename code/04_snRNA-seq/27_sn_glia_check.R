## Louise Huuki-Myers, May 2025
## Compare gilal subtype outputs

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("spatialLIBD")
library("tidyverse")

data_dir <- here("processed-data", "04_snRNA-seq", "27_sn_gila_check")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "27_sn_gila_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

sce <- sce[,sce$cell_type_class == "glia"]
message(Sys.time(), sprintf(" - Subset to %s, ncells = %i", "glia", ncol(sce)))

#### Add cluster info ####

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

ct_cluster_info <- map_dfr(glia_cell_types, function(ct){
    ## load clusters
    cluster_fn <- here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", ct, sprintf("walktrap_snn_k10_subclusters_%s.Rdata", ct))
    load(cluster_fn)
    
    ## load cluster info
    ct_cluster_info <- read_csv(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", ct , sprintf("ERCsn_subtype_clustering_info-%s_k10.csv", ct)),
                                show_col_types = FALSE) |>
        mutate(broad_cell_type = ct)
    
    anno_table <- ct_cluster_info |> select(k10, cell_type_k10) |> column_to_rownames("k10")
    
    sce[,sce$cell_type_broad == ct]$glia_subcluster <<- anno_table[paste0("k", clusters),]

    return(ct_cluster_info)
})

table(sce$glia_subcluster)

table(sce$glia_subcluster, sce$glia_subtype)

## calc concordance
message("Concordance rand score: ", bluster::pairwiseRand(sce$glia_subcluster, sce$glia_subtype, mode = "index"))
# 0.48

#### plot marker genes for subcluster ####
message(Sys.time(), "Plot lit marker genes")
lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    rename(cell_type_target = cell_type) |>
    filter(gene_name %in% rowData(sce)$gene_name)

lit_markers_list <- map(splitit(lit_markers$cell_type_target), ~lit_markers$gene_name[.x])

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_sn_lit_markers_subtype_glia_subcluster.pdf"),
                         cellType_col = "glia_subcluster",
                         gene_name_col = "gene_name")

#### plot QC data ####
pd <- as.data.frame(colData(sce))

qc_cluster_info <- pd |>
                       group_by(!!sym(k_nice), !!sym(cell_type_k)) |>
                       summarize(median_sum = median(sum),
                                 median_detected = median(detected),
                                 median_Mito_percent = median(subsets_Mito_percent),
                                 median_scDblFinder.score = median(scDblFinder.score)) |>
    left_join(cluster_metrics_pass) |>
    arrange(!!sym(cell_type_k))

qc_violin_plot_all <- pd |> 
    left_join(cluster_metrics_pass) |>
    select(!!sym(cell_type_k),  passALL_metricQC, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(!!sym(cell_type_k),  passALL_metricQC), names_to = "metric") |>
    mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))) |>
    ggplot() +
    geom_violin(aes(x = !!sym(cell_type_k), y = value, fill = !!sym(cell_type_k), colour = passALL_metricQC), 
                scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
    # scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    geom_hline(data = cell_class_cutoffs, aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    geom_label(data = cell_class_cutoffs, aes(y = cutoff, label = cutoff_anno), x = 3, color = "blue", vjust = -.5) +
    facet_grid(metric~., scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Quality Metric Value", x = sprintf("subset cell type cluster, %s", k))

ggsave(qc_violin_plot_all, filename = here(plot_dir, sprintf("ERC_sn_QCmetricViolin_subtype_%s-%s.png", cell_type, k_nice)), width = 8, height = 12)


#### sctype ####

## source sc-type functions
message(Sys.time(), "Sctype data")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# prepare gene sets
gs_list <- gene_sets_prepare(db_, "Brain")
rownames(sce) <- rowData(sce)$gene_name

## calc es.max
message(Sys.time(), "- sctype es.max")
es.max <- sctype_score(scRNAseqData = as.matrix(logcounts(sce)),
                       scaled = TRUE,
                       gs = gs_list$gs_positive)

## load previously computed scores
# load(here("processed-data", "04_snRNA-seq", "13_sctype_final", "sctype_es_max-sctype.Rdata"), verbose = TRUE)
## subset to cell type
# es.max <- es.max[,cell_type_index]
dim(es.max)

message(Sys.time(), "- sctype score")
cell_types = cluster_annotation[[cell_type_k]]
# message("cell_types: ", paste(cell_types, collapse = ", "))
cL_results <- purrr::map_dfr(cell_types, function(cluster){
    
    cluster_index = sce[[cell_type_k]] == cluster
    
    es.max.cl = sort(rowSums(es.max[,cluster_index]), decreasing = !0)
    cL_resutls <- head(tibble(cluster = cluster, 
                              sctype = names(es.max.cl), 
                              scores = es.max.cl, 
                              ncells = sum(cluster_index)
    ),10)
    return(cL_resutls)
})

(sctype_score <-  cL_results |> 
    group_by(cluster) |> 
    slice_max(scores)  |>
    mutate(confident = scores >= (ncells/4)) |> 
    mutate(sctype = ifelse(!confident, paste0(sctype, "*"), sctype)))

#### run cell type specific modeling ####
modeling_fn <- here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("modeling_results_subtype-%s_%s.rds", cell_type, k_nice))
pseudobulk_fn = here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("sce_pseudobulk_subtype-%s_%s.rds", cell_type, k_nice))

if(!remodel & file.exists(modeling_fn) & file.exists(pseudobulk_fn)){
    
    message(Sys.time(), " - Load modeling data")
    modeling_results <- readRDS(modeling_fn)
    
} else {
    
    message(Sys.time(), " - Running Spatial Registration on: ", cell_type, " - ", cell_type_k)
    stopifnot(cell_type_k %in% colnames(colData(sce)))
    
    ## Filter to passALL_metricQC nuclei
    sce <- sce[,sce$passALL_metricQC]
    table(sce[[cell_type_k]], sce$sample_id)
    
    ## make APOE syntatic
    sce$APOE <- gsub("/", "", sce$APOE)
    
    modeling_results <-registration_wrapper(
        sce = sce,
        var_registration = cell_type_k,
        var_sample_id = "sample_id",
        covars = c("APOE", "Sex", "Age", "Anc_Afr"),
        gene_ensembl = "gene_id",
        gene_name = "gene_name",
        min_ncells = 10,
        pseudobulk_rds_file = here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("sce_pseudobulk_subtype-%s_%s.rds", cell_type, k_nice))
    )
    
    message(Sys.time(), " - Saving Data")
    saveRDS(modeling_results, file = here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("modeling_results_subtype-%s_%s.rds", cell_type, k_nice)))
    
}

sce_pseudo <- readRDS(here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("sce_pseudobulk_subtype-%s_%s.rds", cell_type, k_nice)))

top_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_genes, here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("subtype_enrichment_top10_%s_%s.csv", cell_type, k_nice)))

#### spatial registration ####

layer_modeling_results <- list()
## sestan sn enrichment
layer_modeling_results$sestan_EC <- readRDS(here("processed-data", "04_snRNA-seq", "24_external_data_check", "sestan_EC_modeling.rds"))
## PsychENCODE enrichment
layer_modeling_results$snDLPFC_PEC <- list()
layer_modeling_results$snDLPFC_PEC$enrichment <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/registration_stats_LIBD.rds")

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
save(spatial_registration, file = here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("ERCsn_subtype_registration-%s_%s.Rdata", cell_type, k_nice)))

# load(here(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check"), sprintf("ERCsn_subtype_registration-%s.Rdata", cell_type)))

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
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_subtype-%s_%s.pdf", ref, cell_type, k_nice)), 
        width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        # reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]]
        # query_colors = cell_type_colors$anno
    ))
    dev.off()
})

#### combine cluster summaries ####
subtype_clustering_info <- 
    cluster_annotation |> 
    left_join(qc_cluster_info) |>
    left_join(sctype_score  |> select(!!sym(cell_type_k) := cluster, sctype)) |>
    left_join(anno_summary  |> rename(!!sym(cell_type_k) := cluster))



# slurmjobs::job_single(name = "26_sn_subtype_check_glia", memory = "50G", command = "Rscript 26_sn_subtype_check.R --cell_type glia", create_shell = TRUE)

# slurmjobs::job_loop(loops = list(cell_type = c("Astro", "Micro", "Endo", "Oligo", "OPC"),
#                                  k = c("10", "20")), 
#                     create_shell = TRUE, 
#                     name = "26_sn_subtype_check", 
#                     create_script = FALSE)

## resubmit Oligo jobs with more mem
# array_submit(
#     "26_sn_subtype_check",
#     task_ids = 6:7,
#     submit = TRUE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



