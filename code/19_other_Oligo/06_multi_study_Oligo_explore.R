## Louise Huuki-Myers, Feb 2026
## Explore multi-study clusters

#### Set up ####
library("SingleCellExperiment")
library("HDF5Array")
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("bluster")
library("ComplexHeatmap")
library("spatialLIBD")

data_dir <- here("processed-data", "19_other_Oligo", "06_multistudy_Oligo_explore")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "06_multistudy_Oligo_explore")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# Import command-line parameters
scec <- matrix(
    c("k", "k", "i", "integer", "k metric for cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

message("load cluster data k=", opt$k)

#### Load multi-study oligo SCE ####
message(Sys.time(), " - Load HDF5 sce") ## loads in 3s
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_build", "sce_multistudy_Oligo"))

## cluster tab
cluster_tab <- readRDS(here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster", sprintf("walktrap_snn_k%02d_subclusters_Multistudy_Oligo.Rds", opt$k)))
# table(cluster_tab$cluster_anno)

sce$Oligo_anno <- cluster_tab$cluster_anno

table(sce[,sce$dataset == "erc"]$Oligo_anno, sce[,sce$dataset == "erc"]$cell_type_fine)

#### Reduced dim plots  - filter to clusters w/ 800+ nuc ####
## examine clsuters w/ over 800 nuc
# length(unique(sce$sample_id)) #85

cluster_anno <- read.csv(here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster", sprintf("Multistudy_Oligo_cluster_annotation_k%i.csv", opt$k)))

cluster_keep <- cluster_anno |> filter(n > 800)
sce <- sce[,sce$Oligo_anno %in% cluster_keep$cluster_anno]
other_oligo_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")


#### compare to original clustering ####
message(Sys.time(), " - Jaccard matrix heatmap")

table(sce$dataset, sce$cell_type_fine)

dataset_jacc_map <- map(c("erc", "dlpfc"), function(d){
    sce_temp <- sce[,sce$dataset == d]
    
    # table(sce_temp$cell_type_fine, sce_temp$Oligo_anno)
    jacc.mat <- linkClustersMatrix(sce_temp$cell_type_fine, sce_temp$Oligo_anno)
    
    ## plot heatmap
    pdf(here(plot_dir, sprintf("Multistudy_Oligo_jaccard_heatmap_k%i_%s.pdf", opt$k, d)), width = 10, height = 10)
    print(Heatmap(jacc.mat,                  
            name = "Correspondence",
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            col = c("black", viridisLite::plasma(100)),
            na_col = "black"))
    dev.off()
    
    return(jacc.mat)
    
})

# #### compare pairwise clustering ####
# k_vals <- c("k10", "k20", "k50")
# names(k_vals) <- k_vals
# cluster_tabs <- map(k_vals, 
#                     ~readRDS(here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster", 
#                                   sprintf("walktrap_snn_%s_subclusters_Multistudy_Oligo.Rds", .x)
#                                   )
#                     )
#                     )

#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("multistudy_Oligo_modeling_results-k%i.rds", opt$k))
pseudobulk_fn = here(data_dir, sprintf("multistudy_Oligo_sce_pseudobulk-k%i.rds", opt$k))

remodel = FALSE

if(!remodel & file.exists(modeling_fn) & file.exists(pseudobulk_fn)){
    
    message(Sys.time(), " - Load modeling data")
    modeling_results <- readRDS(modeling_fn)
    
} else {
    
    message(Sys.time(), " - Running Modeling")
    
    modeling_results <-registration_wrapper(
        sce = sce,
        var_registration = "Oligo_anno",
        var_sample_id = "sample_id",
        covars = c("Age", "Sex"),
        gene_ensembl = "gene_id",
        gene_name = "gene_name",
        min_ncells = 10,
        pseudobulk_rds_file = pseudobulk_fn
    )
    
    message(Sys.time(), " - Saving Data")
    saveRDS(modeling_results, file = modeling_fn)
    
}

sce_pseudo <- readRDS(pseudobulk_fn)

top_enrichment_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_enrichment_genes, here(data_dir, sprintf("multistudy_Oligo_enrichment_top10_k%i.csv", opt$k)), row.names = FALSE)

# #### Register with ERC Oligo subclusters ####
# message(Sys.time(), " - Register vs. ERC Oligo")
# 
# load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)
# 
# erc_oligo_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "modeling_results_subtype-Oligo.rds"))
# 
# cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
#                             modeling_results = erc_oligo_modeling,
#                             model_type = "enrichment",
#                             top_n = 100)
# 
# anno <- annotate_registered_clusters(
#     cor_stats_layer = cor_layer,
#     confidence_threshold = 0.5,
#     cutoff_merge_ratio = 0.1
# )
# 
# ## save data
# cor_layer_anno <- list(cor_layer = cor_layer, anno = anno)
# write_rds(cor_layer_anno, file = here(data_dir, sprintf("multistudy_Oligo_cor_layer_anno_k%i.Rds", opt$k)))
# 
# ## layer stat cor plot
# rownames(cor_layer)[!rownames(cor_layer) %in% names(other_oligo_colors)]
# 
# pdf(here(plot_dir, sprintf("multistudy_Oligo_layer_stat_cor-k%i.pdf", opt$k)))
# layer_stat_cor_plot(cor_layer,
#                     query_colors = other_oligo_colors,
#                     reference_colors = Oligo_OPC_colors,
#                     annotation = anno)
# dev.off()

#### Register with ERC Oligo + OPC subclusters ####
message(Sys.time(), " - Register vs. ERC Oligo + OPC")
if(opt$opc){
    erc_oligoOPC_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "modeling_results_subtype-OligoOPC.rds"))
    
    cor_layer_oligoOPC <- layer_stat_cor(stats = modeling_results$enrichment,
                                         modeling_results = erc_oligoOPC_modeling,
                                         model_type = "enrichment",
                                         top_n = 100)
    
    anno_oligoOPC <- annotate_registered_clusters(
        cor_stats_layer = cor_layer_oligoOPC,
        confidence_threshold = 0.5,
        cutoff_merge_ratio = 0.1
    )
    
    ## save data
    cor_layer_anno <- list(cor_layer = cor_layer_oligoOPC, anno = anno_oligoOPC)
    write_rds(cor_layer_anno, file = here(data_dir, sprintf("multistudy_Oligo_cor_layer_anno_k%i.Rds", opt$k)))
    
    ## layer stat cor plot
    # rownames(cor_layer_oligoOPC)[!rownames(cor_layer_oligoOPC) %in% names(other_oligo_colors)]
    all(rownames(cor_layer_oligoOPC) %in% names(other_oligo_colors))
    
    Oligo_OPC_colors2 <- Oligo_OPC_colors[!grepl("OPC", names(Oligo_OPC_colors))]
    Oligo_OPC_colors2 <- c(Oligo_OPC_colors2, OPC = "#D2B037")
    
    pdf(here(plot_dir, sprintf("multistudy_Oligo_layer_stat_cor-k%i.pdf", opt$k)))
    layer_stat_cor_plot(cor_layer_oligoOPC,
                        query_colors = other_oligo_colors,
                        reference_colors = Oligo_OPC_colors2,
                        annotation = anno_oligoOPC
    )
    dev.off()
}

#### Dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

top_genes_plot <- unique(top_enrichment_genes$gene[top_enrichment_genes$gene %in% rownames(sce)])

all(top_genes_plot %in% rownames(sce))

pdf(here(plot_dir, sprintf("multistudy_Oligo_%s_%s_dotplot_enrichment.pdf", opt$dataset, opt$cluster)), height = 10)
sce |>
    scDotPlot(features = top_genes_plot,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "Marker",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "Marker" = other_oligo_colors),
              scale = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)
sce |>
    scDotPlot(features = top_genes_plot,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "Marker",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "Marker" = other_oligo_colors),
              scale = TRUE,
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

# slurmjobs::job_loop(loops = list(k = c("10", "20", "30")),
#                                  name = "06_multistudy_Oligo_explore",
#                                  create_shell = TRUE,
#                                  create_script = FALSE)

# slurmjobs::array_submit(name = "06_multistudy_Oligo_explore", task_ids = 2, submit = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


