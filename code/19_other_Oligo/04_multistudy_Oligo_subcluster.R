## Louise Huuki-Myers, Jan 2025
## Combine and sub-cluster Oligos from multiple studies 

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

data_dir <- here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_subcluster.R")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "04_multistudy_Oligo_subcluster.R")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

### Load ERC data ####
message(Sys.time(), " - Load ERC sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))
message("ERC n nuc: ", ncol(sce), ", n gene: ", nrow(sce))

sce <- sce[,sce$cell_type_broad %in% c("Oligo", "OPC")]

sce$cell_type_broad <- as.character(sce$cell_type_broad)
table(sce$cell_type_broad)

sce$cell_type_fine <- as.character(sce$cell_type_anno)
table(sce$cell_type_fine)

sce$dataset <- "ERC"

coldata_keep <- c("Barcode", 
                  "sample_id",
                  "dataset",
                  # "chromium_id", 
                  "BrNum", 
                  "Sex", 
                  "Age",
                  # "Diagnosis", 
                  # "Rin",
                  "sum", 
                  "detected", 
                  "subsets_Mito_sum", 
                  "subsets_Mito_detected", 
                  "subsets_Mito_percent", 
                  "total", 
                  # "scDblFinder.sample", 
                  # "scDblFinder.class",
                  "scDblFinder.score",
                  "cell_type_broad",
                  "cell_type_fine")

colData(sce) <- colData(sce)[, coldata_keep]

assay(sce, "binomial_deviance_residuals") <- NULL

#### Load spatialDLPFC ####
message(Sys.time(), " - Load DLPFC sce")

## get DLPFC snRNA-seq data
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce_dlpfc <- HDF5Array::loadHDF5SummarizedExperiment(
    file.path(tempdir(), "sce_DLPFC_annotated")
)
message("DLPFC n nuc: ", ncol(sce_dlpfc), ", n gene:", nrow(sce))

## cell types
sce_dlpfc$cell_type_broad <- as.character(sce_dlpfc$cellType_layer)
sce_dlpfc <- sce_dlpfc[,sce_dlpfc$cell_type_broad %in% c("Oligo", "OPC")]
table(sce_dlpfc$cell_type_broad)

sce_dlpfc$cell_type_fine <- as.character(sce_dlpfc$cellType_hc)
table(sce_dlpfc$cell_type_fine)

## modify colData
# colnames(colData(sce_dlpfc))[colnames(colData(sce_dlpfc)) %in% coldata_keep]

sce_dlpfc$scDblFinder.score <- sce_dlpfc$doubletScore
sce_dlpfc$sample_id <- sce_dlpfc$SAMPLE_ID
sce_dlpfc$Age <- sce_dlpfc$age
sce_dlpfc$Sex <- sce_dlpfc$sex
sce_dlpfc$dataset <- "dlpfc"

all(coldata_keep %in% colnames(colData(sce_dlpfc)))
# coldata_keep[!coldata_keep %in% colnames(colData(sce_dlpfc))]

colData(sce_dlpfc) <- colData(sce_dlpfc)[, coldata_keep]

## modify rowData
rownames(sce_dlpfc) <- rowData(sce_dlpfc)$gene_id

common_genes <- intersect(rownames(sce), rownames(sce_dlpfc))
message("common genes: ", length(common_genes))

sce <- sce[common_genes,]
sce_dlpfc <- sce_dlpfc[common_genes,]
identical(rownames(sce_dlpfc), rownames(sce))

# rowData(sce_dlpfc) <- rowData(sce)

rowRanges(sce_dlpfc) <- NULL 

sce_rd <- rowData(sce)
sce_rr <- rowRanges(sce)

rowRanges(sce) <- NULL 

sce <- cbind(sce, sce_dlpfc)

rowData(sce) <- sce_rd
rownames(sce) <- rownames(sce_rd)

rm(sce_dlpfc)
    
#### Load spatialHPC ####
message(Sys.time(), " - Load HPC sce")

library(ExperimentHub)
ehub <- ExperimentHub()

## Load the HPC dataset
myfiles <- query(ehub, "humanHippocampus2024")
sce_hpc <- myfiles[["EH9606"]]

message("HPC n nuc: ", ncol(sce_hpc), ", n gene:", nrow(sce_hpc))

## cell types
sce_hpc$cell_type_broad <- as.character(sce_hpc$broad.cell.class)
sce_hpc <- sce_hpc[,sce_hpc$cell_type_broad %in% c("Oligo", "OPC")]

## subset to Oligos
sce_hpc$cell_type_fine <- as.character(sce_hpc$superfine.cell.class)

table(sce_hpc$cell_type_fine)
# Oligo.1 Oligo.2 
# 3694    3093

## modify colData
# colnames(colData(sce_hpc))[colnames(colData(sce_hpc)) %in% coldata_keep]

sce_hpc$BrNum <- sce_hpc$brnum
sce_hpc$scDblFinder.score <- sce_hpc$doubletScore
sce_hpc$sample_id <- sce_hpc$Sample
sce_hpc$Age <- sce_hpc$age
sce_hpc$Sex <- sce_hpc$sex
sce_hpc$dataset <- "hpc"

all(coldata_keep %in% colnames(colData(sce_hpc)))
# coldata_keep[!coldata_keep %in% colnames(colData(sce_hpc))]

colData(sce_hpc) <- colData(sce_hpc)[, coldata_keep]

## modify rowData
rownames(sce_hpc) <- rowData(sce_hpc)$gene_id

common_genes <- intersect(rownames(sce), rownames(sce_hpc))
message("common genes: ", length(common_genes))

sce <- sce[common_genes,]
sce_hpc <- sce_hpc[common_genes,]
identical(rownames(sce_hpc), rownames(sce))

rowRanges(sce_hpc) <- NULL 
rowData(sce_hpc) <- rowData(sce)

sce_rd <- rowData(sce)

## match assays
assay(sce, "logcounts") <- NULL
assay(sce, "counts") <- as(assay(sce, "counts"), "Matrix")

assay(sce_hpc, "logcounts") <- NULL
assay(sce_hpc, "counts") <- as(assay(sce_hpc, "counts"), "Matrix")

sce <- cbind(sce, sce_hpc)

rm(sce_hpc)

#### Combined Oligo Data ####

table(sce$dataset)

table(sce$cell_type_broad, sce$dataset)
table(sce$cell_type_fine, sce$dataset)


#### GLM PCA ####
set.seed(425)

harmony_file <- here(data_dir, "multistudy_Oligo_HARMONY_pca.rdata"))

if(file.exists(harmony_file)){
    message(Sys.time(), " - Load saved HARMONY PCs")
    load(harmony_file)
    reducedDim(sce, "HARMONY") <- subtype_HARMONY
} else {
    
    message(Sys.time(), " - running Deviance Feat. Selection")
    sce <- scry::devianceFeatureSelection(sce,
                                          assay = "counts", 
                                          fam = "binomial", 
                                          sorted = F,
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


#### compare vs. previous clusters ####

if(opt$dataset == "spatialDLPFC"){
    
    table(sce$Oligo_anno, sce$cellType_hc)
    jacc.mat <- linkClustersMatrix(sce$Oligo_anno, sce$cellType_hc)
    
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

# slurmjobs::job_single('04_multistudy_Oligo_subcluster', create_shell = TRUE, memory = '100G', command = "Rscript 04_multistudy_Oligo_subcluster.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

