## Louise Huuki-Myers, Jun 2026
## Evaluate label transfer between the snRNA-seq data and the seurat objects for each donor
## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/145f74fb821473fb36d6845c5f6331ea0024f089/code/24_xenium/08_snRNA_xenium_transfer.R

#### Set Up ####

library("sessioninfo")
library("Seurat")
library("here")
library("qs2")

data_dir <- here("processed-data", "25_Seurat", "03_Seurat_label_transfer_evaluation")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("processed-data", "25_Seurat", "03_Seurat_label_transfer_evaluation")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load Data ####
message(Sys.time(), "- Load xenium data")
xen <- qs_read(here("processed-data", "25_Seurat", "02_Seurat_label_transfer", "ERC_seurat_obj_Xenium_label_trasnfer.qs2"))

# xen <- qs_read(here("processed-data", "25_Seurat", "01_Seurat_conversion","ERC_seurat_obj_Xenium.qs2"))

xen

#### Agreement with RCDT ####


#### SEURAT pca ####
head(rownames(xen))
head(rownames(xen))

test <- FeaturePlot(xen, features = c("MOBP"))
