## Louise Huuki-Myers, Feb 2026
## Cluster Multi-study Oligo

#### Set up ####

library("SingleCellExperiment")
library("HDF5Array")
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("scater")
library("scran")
library("igraph")
library("scDotPlot")

data_dir <- here("processed-data", "19_other_Oligo", "05_multistudy_Oligo_cluster")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "05_multistudy_Oligo_cluster")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# Import command-line parameters
scec <- matrix(
    c("k", "k", "i", "integer", "k metric for cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

message("Run clustering with k=", opt$k)

#### Load multi-study oligo SCE ####
message(Sys.time(), " - Load HDF5 sce") ## loads in 3s
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_build", "sce_multistudy_Oligo"))

# message(Sys.time(), " - Load Rds sce") ## 
# sce <- readRDS(here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_build", "sce_multistudy_Oligo.Rds"))

#### SNN + Walktrap cluster ####
set.seed(2526)

## Build SNN graph
message(Sys.time(), " - running buildSNNGraph: k =", opt$k)
snn.gr <- buildSNNGraph(sce, k = opt$k, use.dimred = "HARMONY")

## Run walk trap clustering
message(Sys.time(), " - running walktrap")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters)

saveRDS(clusters, file = here(data_dir, sprintf("Multistudy_Oligo_cluster_k%i.Rds", opt$k)))

# clusters <- readRDS(file = here(data_dir, sprintf("Multistudy_Oligo_cluster_k%i.Rds", opt$k)))

table(clusters)

# clust.50 <- clusterCells(sce, use.dimred="HARMONY", BLUSPARAM=NNGraphParam(k=50))
# table(clust.50)

message(Sys.time(), " - done walktrap")

# #### kmeans clustering #### 
# opt$k = 7
# km <- kmeans(
#     reducedDim(sce, "HARMONY"),
#     centers = opt$k
# )
# table(km$cluster, sce$cell_type_broad)
# table(km$cluster, sce$dataset)
# clusters <- km$cluster

#### Annotate clusters ####
cluster_anno <- as.data.frame(table(clusters, sce$cell_type_broad)) |>
    pivot_wider(names_from = "Var2", values_from = "Freq") |>
    mutate(oligo_prop = Oligo/(Oligo+OPC),
           cell_type = case_when(oligo_prop > 0.8 ~ "Oligo", 
                                 oligo_prop < 0.2 ~ "OPC",
                                 TRUE ~ "Mix"),
           n = Oligo + OPC) |>
    dplyr::rename(cluster = clusters) |>
    group_by(cell_type) |>
    arrange(cell_type, -n) |>
    mutate(cluster_anno = paste0("multi_", cell_type,".", row_number()))

write.csv(cluster_anno, file = here(data_dir, sprintf("Multistudy_Oligo_cluster_annotation_k%i.csv", opt$k)))

cluster_tab <- data.frame(key = colnames(sce), 
                          cluster = clusters, 
                          cluster_anno = cluster_anno$cluster_anno[match(clusters, cluster_anno$cluster)])

head(cluster_tab)

table(cluster_tab$cluster_anno)


#### save cluster_tab data ####
message(Sys.time() , " - saving data")
write_rds(cluster_tab, file = here(data_dir, sprintf("walktrap_snn_k%02d_subclusters_Multistudy_Oligo.Rds", opt$k)))

## Add to sce data & create color pal
sce$Oligo_anno <- cluster_tab$cluster_anno
other_oligo_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")

# sce$metadata$colors <- other_oligo_colors

#### Quality checks ####
pd <- as.data.frame(colData(sce))

colnames(pd)

cell_class_cutoffs <- read_csv(here("processed-data", "04_snRNA-seq", "09_cluster_QC","ERC_sn_cell_class_cutoffs.csv")) |>
    filter(cell_type_class == "glia")

qc_violin_plot_all <- pd |> 
    select(Oligo_anno, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(Oligo_anno), names_to = "metric") |>
    mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))) |>
    ggplot() +
    geom_violin(aes(x = Oligo_anno, y = value, fill = Oligo_anno), 
                draw_quantiles = c(0.25, 0.5, 0.75),
                scale = "width") +
    scale_fill_manual(values = other_oligo_colors)+
    geom_hline(data = cell_class_cutoffs, aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    theme_bw() +
    facet_grid(metric~., scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(qc_violin_plot_all, filename = here(plot_dir, sprintf("Multistudy_Oligo_k%i_QCmetricViolin.png", opt$k)), width = 9)


#### Reduced dim plots ####
message(Sys.time(), " - Plot TSNE")
source(here("code", "utils", "my_plot_reduced_dim.R"))

tsne_plot <- my_plot_reduced_dim(sce,
                                 prefix = "other_Oligo",
                                 dimred = "TSNE",
                                 my_var = "Oligo_anno",
                                 var_type = "cat",
                                 save_plot = TRUE,
                                 suffix = sprintf("multistudy_Oligo_k%i", opt$k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE,
                                 color_pal = other_oligo_colors) +
    theme(legend.position = "None")

ggsave(tsne_plot, filename = here(plot_dir, sprintf("other_Oligo_TSNE-Oligo_anno-multistudy_Oligo_k%i.png", opt$k)))

#### Reduced dim plots  - filter to clusters w/ 800+ nuc ####

## examine clsuters w/ over 800 nuc
# length(unique(sce$sample_id)) #85

cluster_keep <- cluster_anno |> filter(n > 800)

sce <- sce[,sce$Oligo_anno %in% cluster_keep$cluster_anno]

other_oligo_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")

message(Sys.time(), " - Plot TSNE")

tsne_plot <- my_plot_reduced_dim(sce,
                                 prefix = "other_Oligo",
                                 dimred = "TSNE",
                                 my_var = "Oligo_anno",
                                 var_type = "cat",
                                 save_plot = TRUE,
                                 suffix = sprintf("multistudy_Oligo_k%i_filter", opt$k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE,
                                 color_pal = other_oligo_colors) 

# ggsave(tsne_plot, filename = here(plot_dir, sprintf("other_Oligo_TSNE-Oligo_anno-multistudy_Oligo_k%i_filter.png", opt$k)))

#### plot markers ####

rownames(sce) <- rowData(sce)$gene_name

## Oligo lit markers
lit_markers <- list(OPC = c("PDGFRA", "CSPG4", "MAG", "CNP", "A2B5"),
                    Oligo = c("PLP1", "ZFP191", "ZFP488", "ZFP536", "SOX17", "NKX6-2", "SMARCA4", "CD82", "TFR", "MAL"),
                    premyelin_Oligo = c("SOX10", "OLIG1", "OLIG2", "NKX2-2", "CD9", "ENPP6"),
                    myelinating_Oligo = c("BMP4", "ENPP4", "ASAP", "TMEM10", "MOG")
                    # disease_associated = c("SERPINA3", "C4B", "TNFRSF1A", "IL1B", "IL33", "HMOX1", "TNF", "ERK", "ERK2"), #https://doi.org/10.1038/s41593-025-01873-x
                    # AD_risk = c("APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1")
)

lit_markers <- map(lit_markers, ~.x[.x %in% rownames(sce)])
lit_markers <- AnnotationDbi::unlist2(lit_markers)

rowData(sce)$litMarker <- NULL
rowData(sce)$litMarker <- names(lit_markers)[match(rownames(sce), lit_markers)] 
table(rowData(sce)$litMarker)

pdf(here(plot_dir, sprintf("Multistudy_Oligo_subtype_k%i_dotplot_lit.pdf", opt$k)))
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
              scale = FALSE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## ERC Oligo markers

# list.files(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo"))

erc_oligo_enrich <- read.csv(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "subtype_enrichment_top10_OligoOPC.csv"))
# erc_oligo_enrich <- read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_enrichment_top10_Oligo.csv"), row.names = 1)

rowData(sce)$ERC_Oligo <- NULL
rowData(sce)$ERC_Oligo <- erc_oligo_enrich$test[match(rownames(sce), erc_oligo_enrich$gene)] 
table(rowData(sce)$ERC_Oligo)

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

Oligo_OPC_colors <- Oligo_OPC_colors[!grepl("OPC", names(Oligo_OPC_colors))]
Oligo_OPC_colors <- c(Oligo_OPC_colors, OPC = "#D2B037")


pdf(here(plot_dir, sprintf("Multistudy_Oligo_subtype_k%i_dotplot_ERC_Oligo_enrich.pdf", opt$k)), height = 10)
sce |>
    scDotPlot(features = erc_oligo_enrich$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "ERC_Oligo",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "ERC_Oligo" = Oligo_OPC_colors),
              scale = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = erc_oligo_enrich$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "ERC_Oligo",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "ERC_Oligo" = Oligo_OPC_colors),
              scale = TRUE,
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = erc_oligo_enrich$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "ERC_Oligo",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "ERC_Oligo" = Oligo_OPC_colors),
              scale = TRUE,
              clusterRows = TRUE,
              groupLegends = FALSE)

dev.off()


# slurmjobs::job_loop(loops = list(k = c("10", "20", "30")),
#                                  name = "05_multistudy_Oligo_cluster",
#                                  create_shell = TRUE,
#                                  create_script = FALSE)

# slurmjobs::array_submit(name = "05_multistudy_Oligo_cluster", task_ids = 2, submit = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


