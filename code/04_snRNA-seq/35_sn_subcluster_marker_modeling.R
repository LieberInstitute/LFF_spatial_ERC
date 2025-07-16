## Louise Huuki-Myers, January 2025
## Run DeconvoBuddies::get_mean_ratio on snRNA_seq data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("DeconvoBuddies")
library("HDF5Array")
library("tidyverse")
library("getopt")
library("spatialLIBD")

# Import command-line parameters
scec <- matrix(
    c("celltype", "c", "1", "character", "Name of celltype"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

celltype <- opt$celltype

# celltype = "Oligo"

data_dir <- here("processed-data", "04_snRNA-seq", "35_sn_subcluster_markers", celltype)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "35_sn_subcluster_markers", celltype)
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

## subset sce
stopifnot(celltype %in% sce$cell_type_broad)
sce <- sce[,sce$cell_type_broad == celltype]

sce$cell_type_anno <- droplevels(sce$cell_type_anno)
table(sce$cell_type_anno)


#### Run Mean Ratio ####
message(Sys.time(), " - Run MeanRatio on ", celltype)
marker_stats_MeanRatio <- get_mean_ratio(
    sce = sce,
    assay_name = "logcounts",
    cellType_col = "cell_type_anno",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

save(marker_stats_MeanRatio, file = here(data_dir, sprintf("marker_stats_MeanRatio_%s.Rdata", celltype)))
# load(here(data_dir, sprintf("marker_stats_MeanRatio_%s.Rdata", celltype)), verbose = TRUE)

top_MeanRatio_genes <- marker_stats_MeanRatio |> 
    arrange(cellType.target) |>
    filter(MeanRatio > 1, MeanRatio.rank <= 10)

write.csv(top_MeanRatio_genes, here(data_dir, sprintf("subtype_MeanRatio_top10_%s.csv", celltype)), row.names = FALSE)
# top_MeanRatio_genes <- read.csv(here(data_dir, sprintf("subtype_MeanRatio_top10_%s.csv", celltype)))

#### plot top MR markers ####
message(Sys.time(), " - Plots")

cell_type_colors <- metadata(sce)$cell_type_colors

## plot markers
plot_marker_express_ALL(
    sce,
    marker_stats_MeanRatio,
    pdf_fn = here(plot_dir, sprintf("sn_violin_MeanRatio_top10-%s.pdf", celltype)),
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = "cell_type_anno",
    color_pal = cell_type_colors$anno,
    plot_points = FALSE
)

#### MeanRatio dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_MeanRatio_genes$cellType.target[match(rownames(sce), top_MeanRatio_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_MeanRatio.pdf", celltype)))
sce |>
    scDotPlot(features = top_MeanRatio_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "Marker" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()


#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("modeling_results_subtype-%s.rds", celltype))
pseudobulk_fn = here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", celltype))

remodel = FALSE

if(!remodel & file.exists(modeling_fn) & file.exists(pseudobulk_fn)){
    
    message(Sys.time(), " - Load modeling data")
    modeling_results <- readRDS(modeling_fn)
    
} else {
    
    message(Sys.time(), " - Running Modeling")
    ## make APOE syntatic
    sce$APOE <- gsub("/", "", sce$APOE)
    
    modeling_results <-registration_wrapper(
        sce = sce,
        var_registration = "cell_type_anno",
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

sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", celltype)))

top_enrichment_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s.csv", celltype)))


#### Enrichment dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_enrichment.pdf", celltype)))
sce |>
    scDotPlot(features = top_enrichment_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "Marker" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

#### spatial registration ####
# 
# layer_modeling_results <- list()
# 
# ## ERC
# layer_modeling_results$spatialERC <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
# colnames(layer_modeling_results$spatialERC$enrichment) <- gsub("_Sp", "~Sp", colnames(layer_modeling_results$spatialERC$enrichment))
# 
# ## PsychENCODE enrichment
# layer_modeling_results$snDLPFC_PEC <- list()
# layer_modeling_results$snDLPFC_PEC$enrichment <- readRDS("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/registration_stats_LIBD.rds")
# 
# ## sestan sn enrichment
# layer_modeling_results$sestan_EC <- readRDS(here("processed-data", "04_snRNA-seq", "24_external_data_check", "sestan_EC_modeling.rds"))
# 
# ## spatial HPC modeling 
# layer_modeling_results$spatialHPC <- list()
# layer_modeling_results$spatialHPC$enrichment <- read.csv("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)
# 
# names(layer_modeling_results)
# # [1] "HumanPilot"   "spatialDLPFC" "spatialERC"   "snDLPFC_PEC"  "sestan_EC"
# 
# cor_layer <- map(layer_modeling_results, function(layer_mod){
#     cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
#                                 modeling_results = layer_mod,
#                                 model_type = "enrichment",
#                                 top_n = 100)
#     return(cor_layer)
# })
# 
# ## Registration Annotation
# anno <- map(cor_layer, ~annotate_registered_clusters(
#     cor_stats_layer = .x,
#     confidence_threshold = 0.5,
#     cutoff_merge_ratio = 0.1
# ))
# 
# anno_summary <- map2_dfr(anno, names(anno), ~.x |> 
#                              select(cluster, layer_label) |>
#                              mutate(dataset = .y)) |>
#     pivot_wider(names_from = "dataset",
#                 values_from = "layer_label") |>
#     arrange(cluster)
# 
# ## save data
# spatial_registration <- list(cor_layer = cor_layer, anno = anno)
# save(spatial_registration, file = here(data_dir, sprintf("ERCsn_subtype_registration-%s.Rdata", celltype)))
# 
# 
# walk(names(cor_layer), function(ref){
#     message(ref)
#     pdf(here(plot_dir, sprintf("layer_stat_cor_%s_subtype-%s.pdf", ref, celltype)), 
#         width = 6 + (ncol(cor_layer[[ref]])/7))
#     print(layer_stat_cor_plot(
#         cor_stats_layer = cor_layer[[ref]],
#         # reference_colors = layer_colors[[ref]],
#         annotation = anno[[ref]]
#         # query_colors = cell_type_colors$anno
#     ))
#     dev.off()
# })

####  GO analysis ####
library("org.Hs.eg.db")
library("clusterProfiler")

sig_enrichment_genes <- sig_genes_extract(
    n = nrow(modeling_results$enrichment),
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
) |> filter(fdr < 0.001)

sig_enrichment_genes |> count(test)


## ENTREZID look up
entrez_search <- bitr(modeling_results$enrichment$ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez <- sig_enrichment_genes |> 
    left_join(entrez_search, by = c("ensembl" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) 

DE_entrez |> count(test)

## Run GO 
universe <- unique(entrez_search$ENTREZID)
length(universe)

ont <- c("CC","BP","MF")
names(ont) <- ont

go_result <- map(ont, ~compareCluster(ENTREZID ~ test,
                                      data = DE_entrez, 
                                      OrgDb = org.Hs.eg.db,
                                      fun = enrichGO,
                                      universe = universe,
                                      ont = .x, ##ALL,CC,BP,MF
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      # qvalueCutoff = 0.05,
                                      readable = TRUE))

compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))
compare_clus |> count(Cluster)

write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_enrichment_%s.csv", celltype)), row.names = FALSE)

## GO dotplot

pdf(file = here(plot_dir, sprintf("GO_dotplot_enrichment_%s.pdf", celltype)), width = 10, height = 10)
walk2(go_result, names(go_result), 
      ~print(
          dotplot(.x, 
                  x = "Cluster", 
                  showCategory = 3, 
                  label_format = 60)  +
              ggtitle(paste("GO Enrichment:", .y)) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
)
dev.off()


#### Cell type checks ####
## from https://www.biocompare.com/Editorial-Articles/590587-A-Guide-to-Oligodendrocyte-Markers/
oligo_markers <- list(OPC = c("PDGFRA", "CSPG4", "MAG", "CNP", "A2B5"),
                      Oligo = c("PLP1", "ZFP191", "ZFP488", "ZFP536", "SOX17", "NKX6-2", "SMARCA4", "CD82", "TFR", "MAL"),
                      premyelin_Oligo = c("SOX10", "OLIGO1", "OLIGO2", "NKX2-2", "CD9"),
                      myelinating_Oligo = c("BMP4", "ENPP4", "ASAP", "TMEM10", "MOG"))

oligo_markers <- map(oligo_markers, ~.x[.x %in% rownames(sce)])

plot_marker_express_List(
    sce[,sce$cell_type_broad %in% c("Oligo", "OPC")],
    gene_list = oligo_markers,
    cellType_col = "cell_type_anno",
    pdf_fn = here(plot_dir, "sn_violin_Oligo_markers.pdf"),
    color_pal = cell_type_colors$anno
)

## oligo marker dot plot
oligo_markers <- AnnotationDbi::unlist2(oligo_markers)

rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- names(oligo_markers)[match(rownames(sce), oligo_markers)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_lit.pdf", celltype)))
sce |>
    scDotPlot(features = oligo_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()


# slurmjobs::job_single('35_sn_subcluster_markers', create_shell = TRUE, memory = '100G', command = "Rscript 35_sn_subcluster_markers.R -celltype 'ct_broad_k20'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

