## Louise Huuki-Myers, March 2025
## Create Main fig plots from SPE object

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("cowplot")
library("patchwork")
library("DeconvoBuddies")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

plot_dir <- here("plots", "05_spe_correct_cluster", "22_SpD_clean_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

table(spe$SpD, spe$sample_id)

rep_sections_tb <- read.csv(here(data_dir, "rep_section.csv")) |>
    filter(rep_section)

#### load SpD annotations ####

#### Define colors for SpD ####
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

#### Spot plots for representative sections ####

single_vis_clus <- vis_clus(
    spe = spe,
    point_size = 1.7,
    colors = SpD_colors,
    sampleid = "Br5517",
    clustervar = "SpD",
    spatial = FALSE,
    guide_point_size = 5
)

ggsave(single_vis_clus, filename = here(plot_dir, "vis_clus_Br5517.png")) 
ggsave(single_vis_clus, filename = here(plot_dir, "vis_clus_Br5517.pdf")) 

cluster_plots <- map(c("AA", "EA"), function(anc) {
# cluster_plots <- map(unique(spe$APOE), function(anc) {
    
    samples <- rep_sections_tb |> filter(Ancestry == anc) |> arrange(APOE)
    # samples <- rep_sections_tb |> filter(APOE == anc) |> arrange(Ancestry)
    
    cluster_row_plots <- map(samples$sample_id, function(s) {
        vis_clus_plot <- vis_clus(
            spe = spe,
            point_size = 1.5,
            colors = SpD_colors,
            sampleid = s,
            clustervar = "SpD"
        ) +
            labs(title = s)  +
            theme(
                legend.position = "None", ## using heat maps label colors
                axis.title.x = element_blank(),
                text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5)
            )
            
        return(vis_clus_plot)
    })
    cluster_row <- Reduce("+", cluster_row_plots) + plot_layout(nrow = 1)
    # ggsave(cluster_row, filename = here(plot_dir, paste0("vis_clust_",k_label,"_row.png")), width = 18)
    return(cluster_row)
})

cluster_grid <- Reduce("/", cluster_plots)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.pdf"), width = 18, height = 9)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.png"), width = 18, height = 9)


#### Plot reduced dims ####
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(spe,
                          prefix = "ERC_spe",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "SpD",
                          color_pal = SpD_colors))


#### Spatial Registration vs. DLPFC & HPC ####
## get reference layer enrichment statistics
layer_modeling_results <- map(c(HumanPilot = "modeling_results", spatialDLPFC = "spatialDLPFC_Visium_modeling_results"), fetch_data)

## Add spatialDLPFC spatial domain annotations to modeling
dlpfc_anno <- read.csv("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv") |>
    dplyr::filter(bayesSpace == "k09") |>
    mutate(layer_combo2 = gsub(" ", "~", layer_combo2))

dlpfc_colnames <- colnames(layer_modeling_results$spatialDLPFC$enrichment)
pwalk(dlpfc_anno, function(...) dlpfc_colnames <<- gsub(..5, ..3, dlpfc_colnames))
colnames(layer_modeling_results$spatialDLPFC$enrichment) <- dlpfc_colnames
head(layer_modeling_results$spatialDLPFC$enrichment)

## spatial HPC modeling 
layer_modeling_results$spatialHPC <- list()
layer_modeling_results$spatialHPC$enrichment <- read.csv("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)

## ERC Annotated modeling results
erc_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))

colnames(erc_modeling$enrichment) <- gsub("_Sp", "~Sp", colnames(erc_modeling$enrichment))


## correlate 
cor_layer <- map(layer_modeling_results, function(layer_mod){
    
    cor_layer <- layer_stat_cor(stats = erc_modeling$enrichment,
                                modeling_results = layer_mod,
                                model_type = "enrichment",
                                top_n = 100)
    
    cor_layer <- cor_layer[levels(spe$SpD),order(colnames(cor_layer))] ## match factor level for SpD
    
    return(cor_layer)
    
})

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

write_csv(anno_summary, file = here(data_dir, "ERC_SpD_spatial_registration_anno_summary.csv"))


## save data
save(cor_layer, anno, file = here(data_dir, "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"))

## reference colors
layer_colors <- list(HumanPilot = spatialLIBD::libd_layer_colors,
                     spatialDLPFC = NULL) #TODO add spatial domain colors)
                     
map(names(cor_layer), function(ref){
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s.pdf", ref)), width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = SpD_colors,
        cluster_rows = FALSE,
        cluster_columns = FALSE
    ))
    dev.off()
})

map(names(cor_layer), function(ref){
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_cluster.pdf", ref)), width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = SpD_colors
    ))
    dev.off()
})

#### plot top markers ####

## load pseuodbulked data
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds"))

rownames(spe) <- rowData(spe)$gene_name
rownames(spe_pb) <- rowData(spe_pb)$gene_name

top_DEGs <- sig_genes_extract(n = 10,
                              modeling_results = erc_modeling,
                              model_type = "enrichment",
                              sce_layer = spe_pb) |>
    mutate(cellType.target = test)

plot_marker_express_ALL(
    spe_pb,
    top_DEGs,
    pdf_fn = here(plot_dir, "SpD_top_enrichment_genes_pb.pdf"),
    n_genes = 10,
    rank_col = "top",
    anno_col = NULL,
    gene_col = "gene",
    cellType_col = "SpD",
    color_pal = SpD_colors,
    plot_points = TRUE
)

plot_marker_express_ALL(
    spe,
    top_DEGs,
    pdf_fn = here(plot_dir, "SpD_top_enrichment_genes.pdf"),
    n_genes = 10,
    rank_col = "top",
    anno_col = NULL,
    gene_col = "gene",
    cellType_col = "SpD",
    color_pal = SpD_colors,
    plot_points = FALSE
)


# slurmjobs::job_single('22_SpD_clean_plots', create_shell = TRUE, memory = '5G', command = "Rscript 22_SpD_clean_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
