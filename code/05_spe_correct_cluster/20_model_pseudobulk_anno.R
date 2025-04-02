## April 2024, Louise Huuki-Myers
## Run spatialLIBD::registration_wrapper to get layer modeling and pseudobulk data for annotated SpD data

## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("HDF5Array")

data_dir <- here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", cluster_type)
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

cluster_var <- "SpD"

#### Run Spatial Registration Function ####
message(Sys.time(), " - Running Spatial Registration on: ", cluster_var)
stopifnot(cluster_var %in% colnames(colData(spe)))

modeling_results <-registration_wrapper(
    sce = spe,
    var_registration = cluster_var,
    var_sample_id = "sample_id",
    covars = c("APOE", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    # suffix = k_nice,
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, sprintf("spe_pseudobulk-%s.rds", cluster_var))
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("modeling_results-%s.rds", cluster_var)))

#### Extract Top Layer Enrichment Genes ####
sce_pb <- readRDS(here(data_dir, sprintf("spe_pseudobulk-%s.rds", cluster_var)))

top_DEGs <- sig_genes_extract(n = 100,
                              modeling_results = modeling_results,
                              model_type = "enrichment",
                              sce_layer = spe_pb) |>
    left_join(cluster_anno_k9)

write_csv(top_DEGs, file = here(data_dir, "enrichment_modeling_SpD_top100.csv"))


# slurmjobs::job_single('20_model_pseudobulk_anno', create_shell = TRUE, memory = '100G', command = "Rscript 20_model_pseudobulk_anno.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
        
        
