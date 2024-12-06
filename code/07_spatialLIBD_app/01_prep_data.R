## Louise Huuki-Myers, Nov 2024
## Export spe object after removing unneeded images + reduced dims for size and convering assays from HDF5 to sparse Matrix for portability

library("spatialLIBD")
library("HDF5Array")
library("Matrix")
library("here")

## local data
#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

lobstr::obj_size(spe)
# 3.20 GB <- but assays are DelayedArray

#### Convert Assays to sparseMatrix for portability #### 
assayNames(spe)
# [1] "counts"    "logcounts"

message(Sys.time(), " - Convert Assays")
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

#### Extract sig genes #####

modeling_results_k09 <- readRDS(here("processed-data","05_spe_correct_cluster","08_model_pseudobulk","BayesSpace_SVGm", "modeling_results-BayesSpace_SVGm_k09.rds"))
spe_pb_k09 <- readRDS(here("processed-data","05_spe_correct_cluster","08_model_pseudobulk","BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k09.rds"))

##define cluster col - needs better documentation in spatialLIBD::sig_genes_extract_all 

spe_pb_k09$spatialLIBD <- spe_pb_k09$BayesSpace_SVGm_k09

tests <- lapply(modeling_results_k09, function(x) {
    colnames(x)[grep("stat", colnames(x))]
})

k=9
stopifnot(length(tests$anova) == 1) ## assuming only "all"
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))

sig_genes_k09 <- sig_genes_extract_all(
    modeling_results = modeling_results_k09,
    sce_layer = spe_pb_k09,
    n=100
)

saveRDS(sig_genes_k09, file = here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", "sig_genes_k09.rds"))

# slurmjobs::job_single('01_prep_data', create_shell = TRUE, memory = '50G', command = "Rscript 01_prep_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()