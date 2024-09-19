
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
# 5.34 GB

spe
# class: SpatialExperiment 
# dim: 30494 122202 
# metadata(0):
#     assays(2): counts logcounts
# rownames(30494): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(122202): AAACAACGAATAGTTC-1_Br5212 AAACAAGTATCTCCCA-1_Br5212 ... TTGTTTGTATTACACG-1_Br6263
# TTGTTTGTGTAAATTC-1_Br6263
# colData names(76): sample_id in_tissue ... BayesSpace_PCA_Harmony_k27 BayesSpace_PCA_Harmony_k28
# reducedDimNames(5): PCA_p1 PCA_p2 TSNE UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

saveRDS(spe, here("code", "07_spatialLIBD_app", "spe_ERC_app.rds"))
