## Louise Huuki-Myers, April 2025
## Run subtype clustering on sn data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("getopt")
# library("scater")
# library("scran")
# library("scry")

# Import command-line parameters
scec <- matrix(
    c("cell_type", "c", "1", "character", "Name of cell type to sub-cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

cell_type <- opt$cell_type

# cell_type = "Endo"

data_dir <- here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", cell_type)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "26_sn_subtype_check", cell_type)
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

stopifnot(cell_type %in% levels(sce$cell_type_broad))

## subset by Broad cell type
sce <- sce[,sce$cell_type_broad == cell_type]
message(Sys.time(), sprintf(" - Subset to %s, ncells = %i", cell_type, ncol(sce)))

table(sce$sample_id)

## read in cluster data

load(here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", cell_type, sprintf("walktrap_snn_k20_subclusters_%s.Rdata", cell_type)),
     verbose = TRUE) 
# clusters

sce$cell_type_fine <- sprintf("%s.%i", cell_type, clusters)
table(sce$cell_type_fine)

#### run cell type specific modeling ####
## make APOE syntatic
sce$APOE <- gsub("/", "", sce$APOE)

cluster_var <- "cell_type_fine"
message(Sys.time(), " - Running Spatial Registration on: ", cluster_var)
stopifnot(cluster_var %in% colnames(colData(sce)))

table(sce[[cluster_var]], sce$sample_id)

modeling_results <-registration_wrapper(
    sce = sce,
    var_registration = cluster_var,
    var_sample_id = "sample_id",
    covars = c("APOE", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, sprintf("sce_pseudobulk-%s.rds", cluster_var))
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("modeling_results-%s.rds", cluster_var)))


