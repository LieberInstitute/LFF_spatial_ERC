
library("spatialLIBD")
library("HDF5Array")
library("Matrix")
library("here")

## local data
#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

# lobstr::obj_size(spe)
# 3.15 GB <- but assays are DelayedArray

#### Convert Assays to sparaseMatrix for portability #### 
assayNames(spe)
# [1] "counts"    "logcounts"

logcounts(spe) <- as(logcounts(spe), "sparseMatrix")  
counts(spe) <- as(counts(spe), "sparseMatrix")  

## keep 10x UMAP
reducedDimNames(spe)
# [1] "10x_pca"  "10x_tsne" "10x_umap" "PCA_p1"   "PCA_p2"   "TSNE"     "UMAP"     "HARMONY" 

reducedDims(spe)$`10x_pca` <- NULL
reducedDims(spe)$`10x_tsne` <- NULL
reducedDims(spe)$`10x_umap` <- NULL

## Drop images we don't really use in the app
imgData(spe) <- imgData(spe)[
    !imgData(spe)$image_id %in% c("hires", "detected", "aligned"),
]

lobstr::obj_size(spe)
# 2.79 GB

saveRDS(spe, here("code", "07_spatialLIBD_app", "spe_ERC_app.rds"))
