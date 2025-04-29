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
library("DeconvoBuddies")
library(jaffelab)
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

# cell_type = "Astro"

data_dir <- here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", cell_type)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "26_sn_subtype_check", cell_type)
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

stopifnot(cell_type %in% levels(sce$cell_type_broad))

## subset by Broad cell type
cell_type_index <- sce$cell_type_broad == cell_type
sce <- sce[,cell_type_index]
message(Sys.time(), sprintf(" - Subset to %s, ncells = %i", cell_type, ncol(sce)))

## read in cluster data
k_values <- c(k10 = 10, k20 = 20)

clusters <- map_dfc(k_values, ~get(load(here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", cell_type, sprintf("walktrap_snn_k%i_subclusters_%s.Rdata", .x, cell_type)))))

message("number of clusters")
map_int(clusters, max)

message("number clusters n<100")
map_int(clusters, ~sum(table(.x) < 100))

## calc concordance
message("Concordance rand score between subset cluster: ", bluster::pairwiseRand(clusters$k10, clusters$k20, mode = "index"))

message("Concordance rand score between original cluster: ")
map_dbl(clusters, ~bluster::pairwiseRand(sce$cell_type_fine, .x, mode = "index"))


#### Name clusters ####
subtype_k <- paste0("cell_type_k", k_values)
names(subtype_k) <- names(k_values)

# create cell type cluster names ranked on n cells
(cluster_annotations <- map(colnames(clusters), ~clusters |> 
        count(!!sym(.x)) |> 
        arrange(-n) |> 
        mutate(rank = row_number(),
               "cell_type_{.x}" := sprintf("%s.%i", cell_type, rank)))
)

names(cluster_annotations) <- names(subtype_k)

walk(cluster_annotations, function(c){
    clusters <<- clusters |> left_join(c |> select(1, 4))
})

head(clusters)

## add to colData
colData(sce) <- cbind(colData(sce), clusters[, subtype_k])

map(subtype_k, ~table(sce[[.x]]))

#### plot marker genes ####
message(Sys.time(), "Plot lit marker genes")
lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    rename(cell_type_target = cell_type) |>
    filter(gene_name %in% rowData(sce)$gene_name & cell_type_target == cell_type) 

lit_markers_list <- map(splitit(lit_markers$cell_type_target), ~lit_markers$gene_name[.x])

walk(subtype_k, ~plot_marker_express_List(sce, 
                                          lit_markers_list, 
                                          pdf_fn = here(plot_dir, sprintf("ERC_sn_lit_markers_subtype_%s_%s.pdf", cell_type, .x)),
                                          cellType_col = .x,
                                          gene_name_col = "gene_name"))

#### QC check ####

#### sctype ####

## source sc-type functions
message(Sys.time(), "Sctype data")

# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# 
# db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
# 
# message(Sys.time(), "- sctype score")
# es.max <- sctype_score(scRNAseqData = as.matrix(logcounts(sce)), 
#                        scaled = TRUE, 
#                        gs = gs_list$gs_positive)

## load previously computed scores
load(here("processed-data", "04_snRNA-seq", "13_sctype_final", "sctype_es_max-sctype.Rdata"), verbose = TRUE)
## subset to cell type
es.max <- es.max[,cell_type_index]
dim(es.max)

sctype_score <- map2(subtype_k, names(subtype_k), function(cl_col, kname){
    
    cell_types = cluster_annotations[[kname]][[cl_col]]
    # message("cell_types: ", paste(cell_types, collapse = ", "))
    cL_results <- purrr::map_dfr(cell_types, function(cluster){
        
        cluster_index = sce[[cl_col]] == cluster
        
        es.max.cl = sort(rowSums(es.max[,cluster_index]), decreasing = !0)
        cL_resutls <- head(tibble(cluster = cluster, 
                                  type = names(es.max.cl), 
                                  scores = es.max.cl, 
                                  ncells = sum(cluster_index)
        ),10)
        return(cL_resutls)
    })
    
    sctype_scores <-  cL_results |> 
        group_by(cluster) |> 
        slice_max(scores)  |>
        mutate(confident = scores >= (ncells/4))
    
    return(sctype_scores)
})


#### run cell type specific modeling ####
## make APOE syntatic
sce$APOE <- gsub("/", "", sce$APOE)

registration_anno <- map2(subtype_k, names(subtype_k), function(cluster_var, k_name){
    # cluster_var <- "cell_type_k10"
    message(Sys.time(), " - Running Spatial Registration on: ", cell_type, " - ", cluster_var)
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
        pseudobulk_rds_file = here(data_dir, sprintf("sce_pseudobulk_subtype-%s_%s.rds", cell_type, k_name))
    )
    
    message(Sys.time(), " - Saving Data")
    saveRDS(modeling_results, file = here(data_dir, sprintf("modeling_results_subtype-%s_%s.rds", cell_type, k_name)))
    
    sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_subtype-%s_%s.rds", cell_type, k_name)))
    
    top_genes <- sig_genes_extract(
        n = 10,
        modeling_results = modeling_results,
        model_type = "enrichment",
        reverse = FALSE,
        sce_layer = sce_pseudo,
        gene_name = "gene_name"
    )
    
    write.csv(top_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s_%s.csv", cell_type, k_name)))
    
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
    
    ## save data
    spatial_registration <- list(cor_layer = cor_layer, anno = anno)
    save(spatial_registration, file = here(data_dir, sprintf("ERCsn_subtype_registration-%s_%s.Rdata", cell_type, k_name)))
    
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
        pdf(here(plot_dir, sprintf("layer_stat_cor_%s_subtype-%s_%s.pdf", ref, cell_type, k_name)), 
            width = 6 + (ncol(cor_layer[[ref]])/7))
        print(layer_stat_cor_plot(
            cor_stats_layer = cor_layer[[ref]],
            # reference_colors = layer_colors[[ref]],
            annotation = anno[[ref]]
            # query_colors = cell_type_colors$anno
        ))
        dev.off()
    })
    
    return(anno_summary)
})

#### combine cluster summaries ####

cluster_annotations
sctype_score
registration_anno

anno_summary2 <- colData(sce) |> 
    as.data.frame() |> 
    group_by(cluster = cell_type_fine)  |> 
    summarise(n_cells = n(), n_donors = length(unique(sample_id))) |>
    left_join(anno_summary)

write_csv(anno_summary2, file = here(data_dir, sprintf("ERCsn_subtype_registration_anno_summary-%s_%s.csv", cell_type, k_name)))


# slurmjobs::job_loop(loops = list(cell_type = c("Oligo", "Astro", "Micro", "Endo")), create_shell = TRUE, name = "26_sn_subtype_check", create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



