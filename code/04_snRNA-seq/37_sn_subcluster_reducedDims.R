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

# Import command-line parameters
scec <- matrix(
    c("celltype", "c", "1", "character", "Name of celltype"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

celltype <- opt$celltype

# celltype = "Oligo"

data_dir <- here("processed-data", "04_snRNA-seq", "37_sn_subcluster_reducedDims")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "37_sn_subcluster_reducedDims")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### load old data ####

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

sce <- sce[,sce$cell_type_broad == celltype]
any(duplicated(sce$Barcode))

## load Harmony data
harmony_file <- here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", celltype, sprintf("subcluster_HARMONY_pca_%s.rdata", celltype))
file.exists(harmony_file)
load(harmony_file, verbose = TRUE)

# subtype_HARMONY
nrow(subtype_HARMONY) == ncol(sce)
rownames(subtype_HARMONY) <- sce$Barcode


#### Load the updated subcluster data ####
message(Sys.time(), " - Load Updated HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))
sce <- sce[,sce$cell_type_broad == celltype]
any(duplicated(sce$Barcode))

## add Harmony data
rownames(sce) <- rowData(sce)$gene_name
all(sce$Barcode %in% rownames(subtype_HARMONY))

subtype_HARMONY <- subtype_HARMONY[sce$Barcode,]
rownames(subtype_HARMONY) <- colnames(sce)

reducedDims(sce)$HARMONY_st <- subtype_HARMONY

#### clac TSNE & UMAP ####
set.seed(716)

message(Sys.time(), " - run TSNE")
sce <- runTSNE(sce, dimred = "HARMONY_st")

message(Sys.time(), " - run UMAP")
sce <- runUMAP(sce, dimred = "HARMONY_st")

#### plot reduced dims ####
cell_type_colors <- metadata(sce)$cell_type_colors

source(here("code", "utils", "my_plot_reduced_dim.R"))
walk(c("TSNE", "UMAP"),
    ~my_plot_reduced_dim(sce,
                         prefix = "sn_subcluster",
                         dimred = .x,
                         my_var = "cell_type_anno",
                         var_type = "cat",
                         save_plot = TRUE,
                         suffix = celltype,
                         facet = FALSE,
                         plot_dir_rd = plot_dir,
                         verbose = TRUE, 
                         add_label = TRUE)
     )

walk(c("TSNE", "UMAP"),
    ~my_plot_reduced_dim(sce,
                         prefix = "sn_subcluster",
                         dimred = .x,
                         my_var = "cell_type_anno",
                         var_type = "cat",
                         save_plot = TRUE,
                         suffix = paste0(celltype, "_color"),
                         color_pal = cell_type_colors$anno,
                         facet = FALSE,
                         plot_dir_rd = plot_dir,
                         verbose = TRUE, 
                         add_label = TRUE)
     )

# slurmjobs::job_loop(loops = list(celltype = c("Astro", "OPC", "Oligo", "Micro", "Vasc")),
#                     create_shell = TRUE, 
#                     name = "37_sn_subcluster_reducedDims", 
#                     create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




