## Louise Huuki-Myers, Jun 2026
## Convert Xenium and snRNA-seq data to seurat format 
## Adapted from https://github.com/LieberInstitute/spatial_NAc/blob/145f74fb821473fb36d6845c5f6331ea0024f089/code/24_xenium/05_convert_snRNA_seurat.R

#### Set Up ####

library("SingleCellExperiment")
library("sessioninfo")
library("Seurat")
library("here")
library("qs2")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn"
# opt$datatype = "Xenium"

data_dir <- here("processed-data", "25_Seurat", "01_Seurat_conversion")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("processed-data", "25_Seurat", "01_Seurat_conversion")
# if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load experiment  object containing normalized counts ####

if(opt$datatype == "sn"){
    
    message(Sys.time(), " - Load HDF5 sce")
    sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))
    
    message("Remove 0 count genes:", sum(rowSums(assay(sce, "counts")) == 0))
    sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]
    
    message(Sys.time(), "- Convert SCE counts to sparse Matrix")
    counts(sce) <- as(counts(sce), "sparseMatrix")  
    
    message(Sys.time(), "- Convert SCE  logcounts to sparse Matrix")
    logcounts(sce) <- as(logcounts(sce), "sparseMatrix")  
    
    ## set colnames
    colnames(sce) <- sce$Barcode
    
} else if(opt$datatype == "Xenium"){
    
    message(Sys.time(), " - Load QS2 xenium data")
    sce <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
    sce <- sce[,sce$spot_class == "singlet"]
    
    message("Remove 0 count genes:", sum(rowSums(assay(sce, "counts")) == 0))
    sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]
    
    ## set rownames to gene_id
    rownames(sce) <- rowData(sce)$gene_id
    
}

## check for duplicated rownames
message("Any duplicated rownames: ", any(duplicated(rownames(sce))))
message("Any duplicated colnames: ", any(duplicated(colnames(sce))))

# rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
# rownames(sce) <- rowData(sce)$Symbol.uniq

sce

##### Convert to seurat ####
message(Sys.time(), " - Convert to Seurat")

seurat_obj <- as.Seurat(sce,
                        counts = "counts",
                        data = "logcounts")

Assays(seurat_obj)

#### Run PCA ####
message(Sys.time(), " - Run SEURAT PCA")
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj, 
                     seed.use = 1318, 
                     npcs = 50,
                     reduction.name = "seurat_pca")

# print the object
seurat_obj

#### Save the object ####
message(Sys.time(), " - Save qs2")
qs2::qs_save(seurat_obj, file = here(data_dir, sprintf("ERC_seurat_obj_%s.qs2", opt$datatype)))

# slurmjobs::job_single('01_Seurat_conversion_Xenium', create_shell = TRUE, memory = '15G', command = "Rscript 01_Seurat_conversion.R --datatype Xenium")
# slurmjobs::job_single('01_Seurat_conversion_sn', create_shell = TRUE, memory = '50G', command = "Rscript 01_Seurat_conversion.R --datatype sn")


#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()