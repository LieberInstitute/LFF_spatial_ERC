## Louise Huuki-Myers, Dec 2025
## Subcluster Oligos in other datasets, check for Oligo 3

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("scater")
library("scran")
library("scry")
library("DeconvoBuddies")
library("bluster")
library("ComplexHeatmap")
library("getopt")
library("TSCAN")

data_dir <- here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "01_other_Oligo_subcluster")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
scec <- matrix(
    c("dataset", "d", "1", "character", "dataset",
      "k", "k", "2", "integer", "k metric for cluster",
      "ds_short", "ds", "3", "character", "dataset short name",
      "opc", "o", "4", "logical", "dataset short name"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt <- list()

# opt$dataset <- "spatialDLPFC"
# opt$ds_short <- "dlpfc"

# opt$dataset <- "spatialHPC"
# opt$ds_short <- "hpc"

# opt$k <- 10
# opt$opc <- TRUE

message(Sys.time(), " - Data:",  opt$dataset, ", k: ", opt$k, " Include OPC: ", opt$opc)

#### Get Oligo data ####

if(opt$dataset == "spatialDLPFC"){
    
    ## get DLPFC snRNA-seq data
    sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
    sce_path <- unzip(sce_path_zip, exdir = tempdir())
    sce <- HDF5Array::loadHDF5SummarizedExperiment(
        file.path(tempdir(), "sce_DLPFC_annotated")
    )
    
    table(sce$cellType_hc)
    table(sce$cellType_layer == "Oligo")
    table(sce$cellType_layer == "OPC")
    
    sce <- sce[,!is.na(sce$cellType_layer)]
    
    sce$cellType_broad <- sce$cellType_layer
    sce$cellType_fine <- sce$cellType_hc
    
    if(opt$opc){
        ## subset to Oligos & OPC
        sce <- sce[,sce$cellType_layer %in% c("Oligo", "OPC")]
        
        opt$dataset <- paste0(opt$dataset, "_wOPC")
    } else {
        ## subset to just Oligos
        sce <- sce[,sce$cellType_layer == "Oligo"]
    }
    
    sce$cellType_fine <- droplevels(sce$cellType_fine)
    sce$cellType_broad <- droplevels(sce$cellType_broad)
    
    sce$sample_id <- sce$BrNum
    
} else if(opt$dataset == "spatialHPC"){
    
    library(ExperimentHub)
    
    ehub <- ExperimentHub()
    
    ## Load the HPC dataset
    myfiles <- query(ehub, "humanHippocampus2024")
    
    sce <- myfiles[["EH9606"]]
    
    table(sce$broad.cell.class)
    # Neuron  Micro  Astro  Oligo  Other 
    # 52000   3940   4298   6787   8386 
    
    ## subset to Oligos
    sce <- sce[,sce$broad.cell.class == "Oligo"]
    sce$superfine.cell.class <- droplevels(sce$superfine.cell.class)
    
    table(sce$superfine.cell.class)
    # Oligo.1 Oligo.2 
    # 3694    3093
    
    sce$sample_id <- sce$brnum
} else if(opt$dataset == "spatialdACC"){
   load("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/snRNA-seq/05_azimuth/sce_azimuth.Rdata", verbose = TRUE)
    #sce
    colData(sce)
    
    table(sce$cellType_azimuth)
    
    sce$cellType_broad <- sce$cellType_azimuth
    sce$cellType_fine <- sce$cellType_azimuth
    
    if(opt$opc){
        ## subset to Oligos & OPC
        sce <- sce[,sce$cellType_broad %in% c("Oligo", "OPC")]
        
        opt$dataset <- paste0(opt$dataset, "_wOPC")
    } else {
        ## subset to just Oligos
        sce <- sce[,sce$cellType_broad == "Oligo"]
    }
    
}

table(sce$cellType_fine)

#### GLM PCA ####
set.seed(425)

harmony_file <- here(data_dir, sprintf("subcluster_HARMONY_pca_%s.rdata", opt$dataset))

if(file.exists(harmony_file)){
    message(Sys.time(), " - Load saved HARMONY PCs")
    load(harmony_file)
    reducedDim(sce, "HARMONY") <- subtype_HARMONY
} else {
    
    message(Sys.time(), " - running Deviance Feat. Selection")
    sce <- scry::devianceFeatureSelection(sce,
                                          assay = "counts", fam = "binomial", sorted = F,
                                          batch = as.factor(sce$round)
    )
    
    ## plot biononial deviance distribution
    pdf(here(plot_dir, sprintf("sn_reprocess_binomial_deviance_%s.pdf", opt$dataset)))
    plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
         type = "l", xlab = "ranked genes",
         ylab = "binomial deviance", main = "Feature Selection with Deviance"
    )
    abline(v = 2000, lty = 2, col = "red")
    abline(v = 5000, lty = 2, col = "blue")
    dev.off()
    
    message(Sys.time(), " - running nullResiduals")
    sce <- nullResiduals(sce,
                         assay = "counts", 
                         fam = "binomial", # default params
                         type = "deviance"
    )
    
    
    hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:5000]
    
    message(Sys.time(), " - running PCA")
    sce <- runPCA(sce,
                  exprs_values = "binomial_deviance_residuals",
                  subset_row = hdgs,
                  ncomponents = 100,
                  name = "GLMPCA_approx",
                  BSPARAM = BiocSingular::IrlbaParam()
    )
    
    #### Batch Correction ####
    message(Sys.time(), " - HARMONY Correction")
    
    ## Run harmony
    ## needs PCA
    reducedDim(sce, "PCA") <- reducedDim(sce, "GLMPCA_approx")
    
    message("running Harmony - ", Sys.time())
    sce <- harmony::RunHarmony(sce, group.by.vars = "sample_id", verbose = TRUE)
    
    ## Remove redundant PCA
    reducedDim(sce, "PCA") <- NULL
    
    ## save HARMONY dims
    subtype_HARMONY <- reducedDim(sce, "HARMONY")
    save(subtype_HARMONY, file = harmony_file)
}

#### SNN + Walktrap cluster ####

## Build SNN graph
message(Sys.time(), " - running buildSNNGraph: k =", opt$k)
snn.gr <- buildSNNGraph(sce, k = opt$k, use.dimred = "HARMONY")

## Run walk trap clustering
message(Sys.time(), " - running walktrap")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters)

if(opt$opc){
    table(clusters, sce$cellType_broad)
    
    cluster_anno <- as.data.frame(table(clusters, sce$cellType_broad)) |>
        pivot_wider(names_from = "Var2", values_from = "Freq") |>
        mutate(oligo_ratio = Oligo/OPC,
               cell_type = ifelse(oligo_ratio >1, "Oligo", "OPC"),
               n = Oligo + OPC) |>
        dplyr::rename(cluster = clusters) |>
        group_by(cell_type) |>
        arrange(cell_type, -n) |>
        mutate(cluster_anno = paste0(opt$ds_short, "_", cell_type,".", row_number()))
    
} else {
    
    cluster_anno <- stack(table(clusters)) |>
        dplyr::rename(n = values, cluster = ind) |>
        arrange(-n) |>
        mutate(cluster_anno = paste0(opt$ds_short, "_Oligo.", row_number()))
    
}


cluster_tab <- data.frame(key = sce$key, 
                          cluster = clusters, 
                          cluster_anno = cluster_anno$cluster_anno[match(clusters, cluster_anno$cluster)])
head(cluster_tab)

table(cluster_tab$cluster_anno)

## save cluster_tab data
message(Sys.time() , " - saving data")
save(cluster_tab, file = here(data_dir, sprintf("walktrap_snn_k%02d_subclusters_%s.Rdata", opt$k, opt$dataset)))

## add existing cluster annotations
# load(here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster", sprintf("walktrap_snn_k%i_subclusters_%s.Rdata", opt$k, opt$dataset)))
# identical(sce$key, cluster_tab$key)
# sce$Oligo_anno <- cluster_tab$cluster_anno

## Add to sce data & create color pal
sce$Oligo_anno <- cluster_tab$cluster_anno
other_oligo_colors <- create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")

# hpc_Oligo.5 hpc_Oligo.2 hpc_Oligo.4 hpc_Oligo.3 hpc_Oligo.1 
# "#3BB273"   "#FF56AF"   "#663894"   "#F57A00"   "#D2B037" 

#### Quality checks ####
pd <- as.data.frame(colData(sce))

colnames(pd)

qc_violin_plot_all <- pd |> 
    select(Oligo_anno, sum, detected, subsets_Mito_percent, doubletScore)  |>
    pivot_longer(!c(Oligo_anno), names_to = "metric") |>
    mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "doubletScore"))) |>
    ggplot() +
    geom_violin(aes(x = Oligo_anno, y = value, fill = Oligo_anno), 
                draw_quantiles = c(0.25, 0.5, 0.75),
                scale = "width") +
    scale_fill_manual(values = other_oligo_colors) +
    theme_bw() +
    facet_grid(metric~., scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(qc_violin_plot_all, filename = here(plot_dir, sprintf("other_Oligo_%s_k%i_QCmetricViolin.png", opt$dataset, opt$k)))

#### Reduced dim plots ####

message(Sys.time(), " - running TSNE")
sce <- runTSNE(sce, dimred = "HARMONY")
## save
tsne_data <- reducedDim(sce, "TSNE")
write_rds(tsne_data, file = here(data_dir, sprintf("TSNE_k%02d_subclusters_%s.Rds", opt$k, opt$dataset)))

source(here("code", "utils", "my_plot_reduced_dim.R"))

tsne_plot <- my_plot_reduced_dim(sce,
                    prefix = "other_Oligo",
                    dimred = "TSNE",
                    my_var = "Oligo_anno",
                    var_type = "cat",
                    save_plot = TRUE,
                    suffix = sprintf("%s_k%i", opt$dataset, opt$k),
                    facet = FALSE,
                    plot_dir_rd = plot_dir,
                    verbose = TRUE, 
                    add_label = TRUE,
                    color_pal = other_oligo_colors)

tsne_plot_sum_umi <- my_plot_reduced_dim(sce,
                    prefix = "other_Oligo",
                    dimred = "TSNE",
                    my_var = "sum",
                    var_type = "con",
                    save_plot = TRUE,
                    suffix = sprintf("%s_k%i", opt$dataset, opt$k),
                    facet = FALSE,
                    plot_dir_rd = plot_dir,
                    verbose = TRUE, 
                    add_label = TRUE)


#### compare vs. previous clusters ####

    
tsne_plot_og <- my_plot_reduced_dim(sce,
                                    prefix = "other_Oligo",
                                    dimred = "TSNE",
                                    my_var = "cellType_fine",
                                    var_type = "cat",
                                    save_plot = TRUE,
                                    suffix = sprintf("%s_k%i", opt$dataset, opt$k),
                                    facet = FALSE,
                                    plot_dir_rd = plot_dir,
                                    verbose = TRUE, 
                                    add_label = TRUE)


table(sce$Oligo_anno, sce$cellType_hc)
jacc.mat <- linkClustersMatrix(sce$Oligo_anno, sce$cellType_fine)


## plot jaccard matrix
pdf(here(plot_dir, sprintf("other_Oligo_%s_k%i_jaccmat.pdf",opt$dataset, opt$k)))
Heatmap(jacc.mat,
        name = "correspondence",
        col = c("black", viridisLite::plasma(100)),
        na_col = "black",
        column_title = sprintf("%s, k%i", opt$dataset, opt$k)
        )

dev.off()

#### Trajectory analysis ####
message(Sys.time(),  " - Trajectory Anlaysis")
by.cluster <- aggregateAcrossCells(sce, ids=sce$Oligo_anno)
centroids <- reducedDim(by.cluster, "HARMONY")

# Set clusters=NULL as we have already aggregated above.
mst <- createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")

tsne_edge <- my_plot_reduced_dim(sce,
                                 prefix = "other_Oligo",
                                 dimred = "TSNE",
                                 my_var = "Oligo_anno",
                                 var_type = "cat",
                                 save_plot = FALSE,
                                 suffix = sprintf("%s_k%i_traj_edge", opt$dataset, opt$k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE) + 
    geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge), color = "black")

ggsave(tsne_edge, filename = here(plot_dir, sprintf("other_Oligo_TSNE-Oligo_anno_%s_k%i_trajectory.png", opt$dataset, opt$k)))


colLabels(sce) <- sce$Oligo_anno

map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="HARMONY")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

sce$pseudotime <- averagePseudotime(tscan.pseudo) 

tsne_pseudotime <- my_plot_reduced_dim(sce,
                                       prefix = "other_Oligo",
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

ggsave(tsne_pseudotime, filename =here(plot_dir, sprintf("other_Oligo_TSNE-Oligo_anno_%s_k%i_trajectory_pseudotime.png", opt$dataset, opt$k)))

# slurmjobs::job_single('01_other_Oligo_subcluster', create_shell = TRUE, memory = '25G', command = "Rscript 01_other_Oligo_subcluster.R --dataset spatialDLPFC --ds_short dlpfc --k 10 --opc TRUE")
# slurmjobs::job_single('01_other_Oligo_subcluster_dacc', create_shell = TRUE, memory = '25G', command = "Rscript 01_other_Oligo_subcluster.R --dataset spatialdACC --ds_short dacc --k 10 --opc TRUE")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

