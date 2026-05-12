## Louise Huuki-Myers, May 2026
## Run label transfer to label xenium cells with SingleR

#### Set up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("getopt")
library("SingleR")
library("BiocParallel")

# Import command-line parameters
scec <- matrix(
    c("cell_type_col", "c", "1", "character", "cell type column"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
# opt$cell_type_col <- "cell_type_broad"

data_dir <- here("processed-data", "21_Xenium", "09_xenium_label_transfer_singleR")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("plots", "21_Xenium", "08_xenium_label_transfer_RCTD")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## SET UP MULTICORE PARAM
ncores <- 4  # match Sys.getenv('SLURM_CPUS_ON_NODE')
bp <- MulticoreParam(workers = ncores) #or bp <- MulticoreParam(4)

set.seed(202605)

#### load data ####
message(Sys.time(), "- Load data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

## brain logcounts in to mem
logcounts(sce) <- as(logcounts(sce), "dgCMatrix")

cat("NAs:", sum(is.na(assay(sce, "logcounts"))), "\n")
cat("Inf:", sum(is.infinite(assay(sce, "logcounts"))), "\n")
cat("NaN:", sum(is.nan(assay(sce, "logcounts"))), "\n")

rownames(sce) <- rowData(sce)$gene_name

test_genes <- rownames(spe)
ref_genes <- rownames(sce)

cat("Shared genes:", length(intersect(test_genes, ref_genes)), "\n")
cat("Test genes not in ref:", sum(!test_genes %in% ref_genes), "\n")

labels <- sce[[opt$cell_type_col]]
cat("NAs in labels:", sum(is.na(labels)), "\n")
cat("Empty strings:", sum(labels == ""), "\n")
cat("N unique labels:", length(unique(labels)), "\n")

# rows with all zeros can cause issues
all_zero_genes <- assay(sce, "logcounts") == 0
cat("All-zero genes in ref:", sum(all_zero_genes), "\n")
sce <- sce[!all_zero_genes, ]


#### Run SingleR ####
message(Sys.time(), "- Run SingleR")
pred <- SingleR(test = spe[, 1:1000],
                  ref = sce,
                  labels = sce[[opt$cell_type_col]],
                  assay.type.test = "logcounts",
                  num.threads = 1,                  # force single thread in C++ code
                  BPPARAM = BiocParallel::SerialParam()  # force single thread in R
)

# pred <- SingleR(test = spe[,1:1000], 
#                 ref = sce, 
#                 labels = sce[[opt$cell_type_col]], 
#                 assay.type.test="logcounts")

# BSPARAM = BiocSingular::IrlbaParam(),
# BPPARAM = bp

message(Sys.time(), " - Save qs2")
qs2::qs_save(pred, here(data_dir, "SingleR_results_xenium.qs2"))

# slurmjobs::job_single('09_xenium_label_transfer_singleR', create_shell = TRUE, memory = '200G', command = "Rscript 09_xenium_label_transfer_singleR.R --cell_type_col cell_type_anno")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
