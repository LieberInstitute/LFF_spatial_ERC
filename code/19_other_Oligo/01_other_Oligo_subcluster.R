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

data_dir <- here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "01_other_Oligo_subcluster")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
scec <- matrix(
    c("dataset", "d", "1", "character", "dataset",
      "k", "k", "2", "int", "k metric for cluster",
      "ds_short", "ds", "1", "character", "dataset short name"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt <- list()
# opt$dataset <- "spatialHPC"
# opt$k <- 20
# opt$ds_short <- "hpc"

message(Sys.time(), " - Data:",  opt$dataset, ", k: ", opt$k)

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
    
    sce <- sce[,!is.na(sce$cellType_layer)]
    
    ## subset to Oligos
    sce <- sce[,sce$cellType_layer == "Oligo"]
    
    sce$cellType_hc <- droplevels(sce$cellType_hc)
    table(sce$cellType_hc)
    
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
}

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

cluster_anno <- stack(table(clusters)) |>
    dplyr::rename(n = values, cluster = ind) |>
    arrange(-n) |>
    mutate(cluster_anno = paste0(opt$ds_short, "_Oligo.", row_number()))

cluster_tab <- data.frame(key = sce$key, 
                          cluster = clusters, 
                          cluster_anno = cluster_anno$cluster_anno[match(clusters, cluster_anno$cluster)])

head(cluster_tab)

table(cluster_tab$cluster_anno)



## save cluster_tab data
message(Sys.time() , " - saving data")
save(cluster_tab, file = here(data_dir, sprintf("walktrap_snn_k%02d_subclusters_%s.Rdata", opt$k, opt$dataset)))

## Add to sce data & create color pal
sce$Oligo_anno <- cluster_tab$cluster_anno
other_oligo_colors <- create_cell_colors(unique(sce$Oligo_anno))

# hpc_Oligo.5 hpc_Oligo.2 hpc_Oligo.4 hpc_Oligo.3 hpc_Oligo.1 
# "#3BB273"   "#FF56AF"   "#663894"   "#F57A00"   "#D2B037" 

#### Quality checks ####


#### Reduced dim plots ####


#### vs. previous clusters ####

if(opt$dataset == "spatialDLPFC"){
    
    table(sce$Oligo_anno, sce$cellType_hc)
    
    
} else if(opt$dataset == "spatialHPC"){
    
    table(sce$Oligo_anno, sce$superfine.cell.class)
    jacc.mat <- linkClustersMatrix(sce$Oligo_anno, sce$superfine.cell.class)
    
    
    
}

## plot jaccard matrix
pdf(here(plot_dir, sprintf("other_Oligo_%s_k%i_jaccmat.pdf",opt$dataset, opt$k)))
Heatmap(jacc.mat,
        name = "correspondence",
        col = c("black", viridisLite::plasma(100)),
        na_col = "black",
        column_title = sprintf("%s, k%i", opt$dataset, opt$k)
        )

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

