## Louise Huuki-Myers, January 2025
## Run DeconvoBuddies::get_mean_ratio on snRNA_seq data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
# library("DeconvoBuddies")
library("HDF5Array")
library("tidyverse")
library("getopt")
# library("spatialLIBD")
library("scater")
library("TSCAN")

celltype <- "OligoOPC"

data_dir <- here("processed-data", "04_snRNA-seq", "38_sn_subcluster_reducedDims_OligoOPC")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "38_sn_subcluster_reducedDims_OligoOPC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### load old data ####

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

sce <- sce[,sce$celltype_broad %in% c("Oligo", "OPC")]
any(duplicated(sce$Barcode))

#### GLM PCA ####
set.seed(425)

harmony_file <- here(data_dir, sprintf("subcluster_HARMONY_pca_%s.rdata", celltype))

if(file.exists(harmony_file)){
    message(Sys.time(), " - Load saved HARMONY PCs")
    load(harmony_file)
    reducedDim(sce, "HARMONY") <- subtype_HARMONY
} else {
    
    message(Sys.time(), " - running Deviance Feat. Selection")
    sce <- scry::devianceFeatureSelection(sce,
                                          assay = "counts", fam = "binomial", sorted = F,
                                          batch = as.factor(sce$seq_round)
    )
    
    ## plot biononial deviance distribution
    pdf(here(plot_dir, sprintf("sn_reprocess_binomial_deviance_%s.pdf", celltype)))
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


#### clac TSNE & UMAP ####
set.seed(716)

message(Sys.time(), " - run TSNE")
sce <- runTSNE(sce, dimred = "HARMONY")

message(Sys.time(), " - run UMAP")
sce <- runUMAP(sce, dimred = "HARMONY")

#### plot reduced dims ####
celltype_colors <- metadata(sce)$celltype_colors

source(here("code", "utils", "my_plot_reduced_dim.R"))

## plot w/ default colors (easier to tell apart)
walk(c("TSNE", "UMAP"),
     ~my_plot_reduced_dim(sce,
                          prefix = "sn_subcluster",
                          dimred = .x,
                          my_var = "celltype_anno",
                          var_type = "cat",
                          save_plot = TRUE,
                          suffix = celltype,
                          facet = FALSE,
                          plot_dir_rd = plot_dir,
                          verbose = TRUE, 
                          add_label = TRUE)
)

## plot w/ cell type colors
walk(c("TSNE", "UMAP"),
     ~my_plot_reduced_dim(sce,
                          prefix = "sn_subcluster",
                          dimred = .x,
                          my_var = "celltype_anno",
                          var_type = "cat",
                          save_plot = TRUE,
                          suffix = paste0(celltype, "_color"),
                          color_pal = celltype_colors$anno,
                          facet = FALSE,
                          plot_dir_rd = plot_dir,
                          verbose = TRUE, 
                          add_label = TRUE)
)

#### Trajectory analysis ####

by.cluster <- aggregateAcrossCells(sce, ids=sce$celltype_anno)
centroids <- reducedDim(by.cluster, "HARMONY_st")

# Set clusters=NULL as we have already aggregated above.
mst <- createClusterMST(centroids, clusters=NULL)
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")

tsne_edge <- my_plot_reduced_dim(sce,
                                 prefix = "sn_subcluster",
                                 dimred = "TSNE",
                                 my_var = "celltype_anno",
                                 var_type = "cat",
                                 save_plot = FALSE,
                                 suffix = celltype,
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE) + 
    geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge), color = "black")

ggsave(tsne_edge, filename =here(plot_dir, sprintf("sn_subcluster_TSNE-celltype_trajectory-Oligo.png")))


colLabels(sce) <- sce$celltype_anno

map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="HARMONY_st")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

sce$pseudotime <- averagePseudotime(tscan.pseudo) 

tsne_pseudotime <- my_plot_reduced_dim(sce,
                                       prefix = "sn_subcluster",
                                       dimred = "TSNE",
                                       my_var = "pseudotime",
                                       var_type = "con",
                                       save_plot = FALSE,
                                       suffix = celltype,
                                       facet = FALSE,
                                       plot_dir_rd = plot_dir,
                                       verbose = TRUE, 
                                       add_label = TRUE) + 
    geom_line(data=line.data, mapping=aes(x=TSNE1, y=TSNE2, group=edge), color = "black")

ggsave(tsne_pseudotime, filename =here(plot_dir, sprintf("sn_subcluster_TSNE-celltype_pseudotime-Oligo.png")))


# slurmjobs::job_single('38_sn_subcluster_reducedDims_OligoOPC', create_shell = TRUE, memory = '25G', command = "Rscript 38_sn_subcluster_reducedDims_OligoOPC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




