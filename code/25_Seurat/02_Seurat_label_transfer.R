## Louise Huuki-Myers, Jun 2026
## Perform label transfer between the snRNA-seq data and the seurat objects for each donor
## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/145f74fb821473fb36d6845c5f6331ea0024f089/code/24_xenium/08_snRNA_xenium_transfer.R

#### Set Up ####

library("sessioninfo")
library("Seurat")
library("here")
library("qs2")


data_dir <- here("processed-data", "25_Seurat", "02_Seurat_label_transfer")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("processed-data", "25_Seurat", "02_Seurat_label_transfer")
# if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load Data ####
message(Sys.time(), "- Load xenium data")
xen <- qs_read(here("processed-data", "25_Seurat", "01_Seurat_conversion","ERC_seurat_obj_Xenium.qs2"))

message(Sys.time(), "- Load snRNA data")
snRNA <- qs_read(here("processed-data", "25_Seurat", "01_Seurat_conversion","ERC_seurat_obj_Xenium.qs2"))

#### Perform the label transfer ####

message(Sys.time(), "- Find transfer anchors")
anchors <- FindTransferAnchors(reference = snRNA, 
                               query = xen, 
                               dims = 1:30,
                               reference.reduction = "seurat_pca",
                               reduction = "rpca")

message(Sys.time(), "- Label Transfer")
predictions <- TransferData(anchorset = anchors, 
                            refdata = snRNA$cell_type_anno, 
                            dims = 1:30)

#Add predictions to the object
message(Sys.time(), "Adding metadata")
xen <- AddMetaData(xen, 
                   metadata = predictions)

#### Save ####

#Save the predictions  
message(paste0("Saving objects - ", Sys.time()))
saveRDS(object = predictions, file = here(data_dir, "xenium_snRNA_predictions.Rds"))

## save updated xen object
qs2::qs_save(xen, file = here(data_dir, "ERC_seurat_obj_Xenium_label_trasnfer.qs2"))

#### Reproduciblity ####

slurmjobs::job_single('02_Seurat_label_transfer', create_shell = TRUE, memory = '15G', command = "Rscript 02_Seurat_label_transfer.R")


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()