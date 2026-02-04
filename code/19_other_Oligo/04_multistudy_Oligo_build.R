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
library("ExperimentHub")
library("HDF5Array")

data_dir <- here("processed-data", "19_other_Oligo", "04_multistudy_Oligo_build")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "04_multistudy_Oligo_build")
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

sce$dataset <- "erc"

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

## add dataset to barcode
colnames(sce) <- paste0(sce$dataset, "_", sce$Barcode)
head(colnames(sce))

## start base for combining data

all_counts <- assay(sce, "counts")
all_counts[1:5, 1:5]
# colnames(all_counts) <- colnames(sce)

all_colData <- colData(sce)

all_rowData <- rowData(sce)

message("Start All data")
message("n nuc matches: ", identical(ncol(all_counts), nrow(all_colData)))
message("colnames match: ", identical(colnames(all_counts), rownames(all_colData)))
message("rownames match: ", identical(rownames(all_counts), rownames(all_rowData)))

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

br_n10 <- as.data.frame(unique(colData(sce_dlpfc)[,c("BrNum","Age", "Sex")]))
rownames(br_n10) <- br_n10$BrNum

all(coldata_keep %in% colnames(colData(sce_dlpfc)))
# coldata_keep[!coldata_keep %in% colnames(colData(sce_dlpfc))]

## add dataset to barcode
colnames(sce_dlpfc) <- paste0(sce_dlpfc$dataset, "_", colnames(sce_dlpfc))
head(colnames(sce_dlpfc))

## modify rowData
rownames(sce_dlpfc) <- rowData(sce_dlpfc)$gene_id

common_genes <- intersect(rownames(all_counts), rownames(sce_dlpfc))
message("common genes: ", length(common_genes))

## combine data
message("merge DLPFC data")
all_counts <- cbind(all_counts[common_genes,], assay(sce_dlpfc, "counts")[common_genes,])
all_colData <- rbind(all_colData, colData(sce_dlpfc)[, coldata_keep])


message("n nuc matches: ", identical(ncol(all_counts), nrow(all_colData)))
message("colnames match: ", identical(colnames(all_counts), rownames(all_colData)))

rm(sce_dlpfc)

#### Load spatialHPC ####
message(Sys.time(), " - Load HPC sce")

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

## add dataset to barcode
colnames(sce_hpc) <- paste0(sce_hpc$dataset, "_", colnames(sce_hpc))
head(colnames(sce_hpc))


## modify rowData
rownames(sce_hpc) <- rowData(sce_hpc)$gene_id

common_genes <- intersect(rownames(all_counts), rownames(sce_hpc))
message("common genes: ", length(common_genes))

## combine data
message("merge HPC data")
all_counts <- cbind(as(all_counts[common_genes,], "Matrix"), 
                    assay(sce_hpc, "counts")[common_genes,])

all_colData <- rbind(all_colData, colData(sce_hpc)[, coldata_keep])

message("n nuc matches: ", identical(ncol(all_counts), nrow(all_colData)))
message("colnames match: ", identical(colnames(all_counts), rownames(all_colData)))

rm(sce_hpc)

#### Load dACC data ####
load("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/snRNA-seq/05_azimuth/sce_azimuth.Rdata", verbose = TRUE)
message("dacc n nuc: ", ncol(sce), ", n gene:", nrow(sce))

## cell types
sce$cell_type_broad <- as.character(sce$cellType_azimuth)
sce <- sce[,sce$cell_type_broad %in% c("Oligo", "OPC")]

## subset to Oligos
sce$cell_type_fine <- as.character(sce$cellType_azimuth)

table(sce$cell_type_fine)
# Oligo.1 Oligo.2 
# 3694    3093

## modify colData
# colnames(colData(sce))[colnames(colData(sce)) %in% coldata_keep]

sce$BrNum <- sce$brain
sce$sample_id <- sce$Sample

all(sce$BrNum %in% br_n10$BrNum)

sce$Age <- br_n10[sce$BrNum,"Age"]
sce$Sex <- br_n10[sce$BrNum,"Sex"]
sce$dataset <- "dacc"

all(coldata_keep %in% colnames(colData(sce)))

## add dataset to barcode
colnames(sce) <- paste0(sce$dataset, "_", colnames(sce))
head(colnames(sce))

## modify rowData
rownames(sce) <- rowData(sce)$gene_id

common_genes <- intersect(rownames(all_counts), rownames(sce))
message("common genes: ", length(common_genes))

## combine data
message("merge dacc data")
all_counts <- cbind(all_counts[common_genes,],
                    assay(sce, "counts")[common_genes,])

all_colData <- rbind(all_colData, colData(sce)[, coldata_keep])

message("n nuc matches: ", identical(ncol(all_counts), nrow(all_colData)))
message("colnames match: ", identical(colnames(all_counts), rownames(all_colData)))

sce$key <- colnames(sce)

#### Combined Oligo Data ####
message(Sys.time(), " - Build Multi-Study Oligo SCE")
sce <- SingleCellExperiment(assay = list(counts = all_counts),
                            colData = all_colData,
                            rowData = all_rowData[common_genes,])

message("ALL data n nuc: ", ncol(sce), ", n gene:", nrow(sce))

table(sce$dataset)

table(sce$cell_type_broad, sce$dataset)
table(sce$cell_type_fine, sce$dataset)

sce$cell_type_dataset <- paste0(sce$cell_type_broad, "_",sce$dataset)

table(sce$cell_type_dataset)

## rm unneeded data
rm(all_counts)

#### GLM PCA ####
set.seed(425)

harmony_file <- here(data_dir, "multistudy_Oligo_HARMONY_pca.rdata")


message(Sys.time(), " - running Deviance Feat. Selection")
sce <- scry::devianceFeatureSelection(sce,
                                      assay = "counts", 
                                      fam = "binomial", 
                                      sorted = F,
                                      batch = as.factor(sce$dataset)
)

## plot biononial deviance distribution
pdf(here(plot_dir, "multistudy_Oligo_binomial_deviance.pdf"))
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

## get logCounts
sce <- logNormCounts(sce)

message(Sys.time(), " - running PCA")
sce <- scater::runPCA(sce,
                      exprs_values = "binomial_deviance_residuals",
                      subset_row = hdgs,
                      ncomponents = 100,
                      name = "GLMPCA_approx",
                      BSPARAM = BiocSingular::IrlbaParam()
)

reducedDimNames(sce)

## un-corrected TSNE 
message(Sys.time(), " - running TSNE")
sce <- runTSNE(sce, dimred = "GLMPCA_approx")

source(here("code", "utils", "my_plot_reduced_dim.R"))

tsne_plot <- my_plot_reduced_dim(sce,
                                 prefix = "Multistudy_Oligo",
                                 dimred = "TSNE",
                                 my_var = "cell_type_dataset",
                                 var_type = "cat",
                                 save_plot = TRUE,
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE)

TSNE_uncorrected <- reducedDims(sce, "TSNE")
saveRDS(TSNE_uncorrected, file = here(data_dir, "Multistudy_Oligo_TSNE_uncorrected.Rds"))

#### Batch Correction ####
message(Sys.time(), " - HARMONY Correction")

## Run harmony
## needs PCA
reducedDim(sce, "PCA") <- reducedDim(sce, "GLMPCA_approx")

message("running Harmony - ", Sys.time())
sce <- harmony::RunHarmony(sce, group.by.vars = "dataset", verbose = TRUE)

## Remove redundant PCA
reducedDim(sce, "PCA") <- NULL

## save HARMONY dims
subtype_HARMONY <- reducedDim(sce, "HARMONY")
save(subtype_HARMONY, file = harmony_file)


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
                                 suffix = sprintf("mulristudy_Oligo_k%i", k),
                                 facet = FALSE,
                                 plot_dir_rd = plot_dir,
                                 verbose = TRUE, 
                                 add_label = TRUE,
                                 color_pal = other_oligo_colors)

#### Save data ####
message(Sys.time(), " - Save HDF5 Data")
saveHDF5SummarizedExperiment(sce, dir = here(data_dir, "sce_multistudy_Oligo"), replace=TRUE)


# slurmjobs::job_single('04_multistudy_Oligo_build', create_shell = TRUE, memory = '100G', command = "Rscript 04_multistudy_Oligo_build.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

