## Louise Huuki-Myers, April 2025
## Pseudobulk single nuc dataset for DGE

library("spatialLIBD")
library("SingleCellExperiment")
library("HDF5Array")
library("dplyr")
library("here")
library("sessioninfo")
library("getopt")
library("scater")

# Import command-line parameters
scec <- matrix(
    c(  "cluster", "c", "1", "character", "Name of cluster"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

## for test
# opt$cluster <- "cell_type_anno"

#### Set up dirs ####
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))
stopifnot(opt$cluster  %in% colnames(colData(sce)))

# rowData(sce)

## add syntactic APOE vars
sce$APOE_syn <- factor(gsub("/", ".", sce$APOE))
levels(sce$APOE_syn)

sce$APOE_carrier_syn <- factor(gsub("\\+", "", sce$APOE_carrier))
levels(sce$APOE_carrier_syn)

## add continuous APOE variable
# APOE_num <- data.frame(APOE_num = 0:3)
# rownames(APOE_num) <- c("E2/E2", "E2/E3", "E3/E4", "E4/E4")
# sce$APOE_num <- APOE_num[sce$APOE, ] 
# table(sce$APOE_num, sce$APOE)

## get mito genes
mito_genes <- as.logical(seqnames(sce) == "chrM")
table(mito_genes)

table(sce[[opt$cluster]])

#### run pseudobulk - cell_type_anno ####
message(Sys.time(), "pseudobulk ", opt$cluster)
sce_pseudo <- registration_pseudobulk(
    sce,
    var_registration = opt$cluster,
    var_sample_id = "sample_id",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL,
    filter_expr = FALSE,
    mito_gene = mito_genes
)

message(Sys.time(), " - Done pseudobulk")

message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))

#### Check n samples for each cell type ####
table(sce_pseudo$APOE_carrier, sce_pseudo$cell_type_broad)
table(sce_pseudo$APOE_carrier, sce_pseudo$registration_variable)

cell_type_count <- colData(sce_pseudo) |> as.data.frame() |> dplyr::count(APOE_carrier, registration_variable)

## cell type must have two or more samples on either side of DEG split (APOE carrier)
enough_samples <- cell_type_count |>
    filter(n >=2) |>
    dplyr::count(registration_variable) |>
    filter(n >=2) |>
    dplyr::pull(registration_variable)

message("Too few samples in: ", 
        paste(levels(sce_pseudo$registration_variable)[!levels(sce_pseudo$registration_variable) %in% enough_samples], collapse = ", ")
)

cell_type_count |>
    mutate(enough_samples = registration_variable %in% enough_samples)|>
    write.csv(here(data_dir, sprintf("sn_psuedobulk_sample_count-%s.csv", opt$cluster)))

## drop too few sample cell types
sce_pseudo <- sce_pseudo[, sce_pseudo$registration_variable %in% enough_samples]
sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)

#### Add PCAs ####
sce_pseudo <- scater::runPCA(sce_pseudo, 
                             ncomponents = 50,
                             name = "PCA")


#### Additional edits + Save ####
## drop all NA cols
all_na <- sapply(colData(sce_pseudo), function(x)all(is.na(x)))
colData(sce_pseudo) <- colData(sce_pseudo)[, names(all_na)[!all_na]]

## save 
message(Sys.time(), " - Save")
saveRDS(sce_pseudo, file = here(data_dir, sprintf("sce_pseudo_DGE-%s.RDS", opt$cluster)))

# sce_pseudo <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn","01_pseudobulk_data", "sce_pseudo_DGE.RDS"))

# slurmjobs::job_single('01_pseudobulk_data_sn', create_shell = TRUE, memory = '100G', command = "Rscript 01_pseudobulk_data.R")

# slurmjobs::job_loop(loops = list(cluster = c("cell_type_anno", "cell_type_broad")),
#                     name = "01_pseudobulk_data_sn",
#                     create_shell = TRUE,
#                     create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
