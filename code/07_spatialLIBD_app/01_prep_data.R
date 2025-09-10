## Louise Huuki-Myers, Nov 2024
## Export spe object after removing unneeded images + reduced dims for size and convering assays from HDF5 to sparse Matrix for portability

library("spatialLIBD")
library("HDF5Array")
library("Matrix")
library("here")
library("sessioninfo")

## local data
#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

## spatialLIBD col
spe$spatialLIBD <- spe$SpD

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

modeling_results <- readRDS(here("processed-data","05_spe_correct_cluster","20_model_pseudobulk_anno","modeling_results-SpD.rds"))
spe_pb <- readRDS(here("processed-data","05_spe_correct_cluster","20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds"))

colnames(rowData(spe_pb))

##define cluster col - needs better documentation in spatialLIBD::sig_genes_extract_all 
spe_pb$spatialLIBD <- spe_pb$SpD

## Add gene_search
# rowData(spe_pb)$gene_search <- paste0(rowData(spe_pb)$gene_name, "; ", rowData(spe_pb)$gene_id)

tests <- lapply(modeling_results, function(x) {
    colnames(x)[grep("stat", colnames(x))]
})

# names(tests)

# modeling_results$anova

k=9
stopifnot(length(tests$anova) == 1) ## assuming only "all"
stopifnot(length(tests$enrichment) == k)
stopifnot(length(tests$pairwise) == choose(k, 2))


# check_sce_layer(spe_pb)
# check_modeling_results(modeling_results)

sig_genes <- sig_genes_extract_all(
    modeling_results = modeling_results,
    sce_layer = spe_pb,
    n = nrow(spe_pb),
    gene_name = "gene_id"
)

nrow(sig_genes)
# 1023032

saveRDS(sig_genes, file = here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "sig_genes_SpD.rds"))

# slurmjobs::job_single('01_prep_data', create_shell = TRUE, memory = '50G', command = "Rscript 01_prep_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()