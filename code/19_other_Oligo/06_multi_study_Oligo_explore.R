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
library("scDotPlot")
library("scater")
library("scran")
library("TSCAN")

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
rownames(sce) <- rowData(sce)$gene_name

## cluster tab
cluster_tab <- readRDS(here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster", sprintf("walktrap_snn_k%02d_subclusters_Multistudy_Oligo.Rds", opt$k)))
# table(cluster_tab$cluster_anno)

sce$Oligo_anno <- cluster_tab$cluster_anno
sce$Oligo_broad <- gsub("multi_|.[0-9]+","", sce$Oligo_anno)

table(sce$Oligo_broad, sce$cell_type_broad)
    
table(sce[,sce$dataset == "erc"]$Oligo_anno, sce[,sce$dataset == "erc"]$cell_type_fine)

#### Reduced dim plots  - filter to clusters w/ 800+ nuc ####
## examine clsuters w/ over 800 nuc
# length(unique(sce$sample_id)) #85

cluster_anno <- read.csv(here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster", sprintf("Multistudy_Oligo_cluster_annotation_k%i.csv", opt$k)))

cluster_keep <- cluster_anno |> filter(n > 800)
sce <- sce[,sce$Oligo_anno %in% cluster_keep$cluster_anno]
other_oligo_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")

table(sce[,sce$dataset == "erc"]$Oligo_anno, sce[,sce$dataset == "erc"]$cell_type_fine)

#### proportions ####

oligo_type_prop <- colData(sce) |>
    as.data.frame() |>
    group_by(dataset, Oligo_broad, Oligo_anno) |>
    summarize(n = n()) |>
    group_by(dataset, Oligo_broad) |>
    mutate(prop = n/sum(n))

n_nuclei_barplot <- oligo_type_prop |>
    ggplot(aes(x = dataset, y = n, fill = Oligo_anno)) +
    geom_col() + 
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2) +
    facet_wrap(~Oligo_broad) +
    scale_fill_manual(values = other_oligo_colors) +
    theme_bw()

ggsave(n_nuclei_barplot, filename = here(plot_dir, "multistudy_Oligo_n_barplot.png"))

prop_barplot <- oligo_type_prop |>
    ggplot(aes(x = dataset, y = prop, fill = Oligo_anno)) +
    geom_col() + 
    geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5), size = 2) +
    facet_wrap(~Oligo_broad) +
    scale_fill_manual(values = other_oligo_colors)+
    theme_bw()

ggsave(prop_barplot, filename = here(plot_dir, "multistudy_Oligo_prop_barplot.png"))


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

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

Oligo_OPC_colors2 <- Oligo_OPC_colors[!grepl("OPC", names(Oligo_OPC_colors))]
Oligo_OPC_colors2 <- c(Oligo_OPC_colors2, OPC = "#D2B037")

pdf(here(plot_dir, sprintf("multistudy_Oligo_layer_stat_cor-k%i.pdf", opt$k)))
layer_stat_cor_plot(cor_layer_oligoOPC,
                    query_colors = other_oligo_colors,
                    reference_colors = Oligo_OPC_colors2,
                    annotation = anno_oligoOPC
)
dev.off()

#### Dot plots ####


## Oligo lit markers 2.0
lit_markers <- consensus_marker_sets <- list(
    OPC = c("PDGFRA", "PTPRZ1", "CSPG4", "VCAN", "CNR1", "COL14A1"),
    NFOL = c("GPR17", "BMP4", "CD9", "TMEM2", "ARPC1B"),
    MOL = c("OPALIN", "MOBP", "MAL", "PLP1", "APOD", "SEPP1", "S100B"),
    Reactive_OL = c("C4B", "SERPINA3N", "CD74", "STAT1", "IRF7","IL33", "BST2", "IFIT1")
)

lit_markers <- map(lit_markers, ~.x[.x %in% rownames(sce)])
lit_markers <- AnnotationDbi::unlist2(lit_markers)

rowData(sce)$litMarker <- NULL
rowData(sce)$litMarker <- names(lit_markers)[match(rownames(sce), lit_markers)] 
table(rowData(sce)$litMarker)

pdf(here(plot_dir, sprintf("Multistudy_Oligo_subtype_k%i_dotplot_lit2.pdf", opt$k)))
sce |>
    scDotPlot(features = lit_markers,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "litMarker",
              scale = TRUE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = lit_markers,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "litMarker",
              scale = TRUE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = TRUE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = lit_markers,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "litMarker",
              scale = FALSE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)

dev.off()


rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

top_genes_plot <- unique(top_enrichment_genes$gene[top_enrichment_genes$gene %in% rownames(sce)])

all(top_genes_plot %in% rownames(sce))

pdf(here(plot_dir, sprintf("multistudy_Oligo_k%i_dotplot_enrichment.pdf", opt$k)), height = 12)
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
sce |>
    scDotPlot(features = top_genes_plot,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "Marker",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "Marker" = other_oligo_colors),
              scale = TRUE,
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()

#### Trajectory analysis ####
message(Sys.time(),  " - Trajectory Anlaysis")
by.cluster <- scater::aggregateAcrossCells(sce, ids=sce$Oligo_anno)
centroids <- reducedDim(by.cluster, "HARMONY")

# Set clusters=NULL as we have already aggregated above.
mst <- TSCAN::createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")
source(here("code", "utils", "my_plot_reduced_dim.R"))

tsne_edge <- my_plot_reduced_dim(sce,
                                 prefix = "multistudy_Oligo",
                                 dimred = "TSNE",
                                 my_var = "Oligo_anno",
                                 var_type = "cat",
                                 save_plot = FALSE,
                                 suffix = sprintf("k%i_traj_edge", opt$k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE) + 
    geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge), color = "black")

ggsave(tsne_edge, filename = here(plot_dir, sprintf("multistudy_Oligo_TSNE-Oligo_anno_k%i_trajectory.png", opt$k)))


colLabels(sce) <- sce$Oligo_anno

map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="HARMONY")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

sce$pseudotime <- averagePseudotime(tscan.pseudo) 

tsne_pseudotime <- my_plot_reduced_dim(sce,
                                       prefix = "multistudy_Oligo",
                                       dimred = "TSNE",
                                       my_var = "pseudotime",
                                       var_type = "con",
                                       save_plot = FALSE,
                                       suffix = data,
                                       facet = FALSE,
                                       plot_dir_rd = plot_dir,
                                       verbose = TRUE, 
                                       add_label = TRUE) + 
    geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge), color = "black")

ggsave(tsne_pseudotime, filename =here(plot_dir, sprintf("multistudy_Oligo_TSNE-Oligo_anno_k%i_trajectory_pseudotime.png", opt$k)))

#### Hierarchial clustering ####
by.cluster <- logNormCounts(by.cluster)
# Error in .local(x, ...) : size factors should be positive

## Perform hierarchical clustering
message(Sys.time(), " - Cluster Again")
dist.clusCollapsed <- dist(t(counts(by.cluster)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)

# Print for future reference
pdf(here(plot_dir, "multistudy_Oligo_dend.pdf"), height = 4)
par(cex = 0.6, font = 2)

## cell type colors on leaves & branches
dend |> 
    # set("leaves_col", cell_type_colors$anno[labels(dend)]) |> 
    # set("leaves_cex", 2) |> 
    # set("leaves_pch", 19) |> 
    # set("branches_k_color", k = 7,
    #     value = cell_type_colors$broad[c("Inhib", "Excit", "Micro", "Vasc", "Oligo", "Astro", "OPC")]) |>
    plot()

# abline(v = 500, lty = 2)
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


