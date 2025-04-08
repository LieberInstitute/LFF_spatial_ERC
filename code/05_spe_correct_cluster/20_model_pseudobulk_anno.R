## April 2024, Louise Huuki-Myers
## Run spatialLIBD::registration_wrapper to get layer modeling and pseudobulk data for annotated SpD data

## Required libraries
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("HDF5Array")

data_dir <- here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

## add syntacticly valid version of SpD
spe$SpD_syn <- gsub("~", "_", spe$SpD)
table(spe$SpD_syn)

## make APOE syntatic
spe$APOE <- gsub("/", "", spe$APOE)

#### Run Spatial Registration Function ####
message(Sys.time(), " - Running Spatial Registration on: SpD_syn")

modeling_results <-registration_wrapper(
    sce = spe,
    var_registration = "SpD_syn",
    var_sample_id = "sample_id",
    covars = c("APOE", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    # suffix = k_nice,
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, "spe_pseudobulk-SpD.rds")
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, "modeling_results-SpD.rds"))

#### Extract Top Layer Enrichment Genes ####
spe_pb <- readRDS(here(data_dir, "spe_pseudobulk-SpD.rds"))

top_DEGs <- sig_genes_extract(n = 100,
                              modeling_results = modeling_results,
                              model_type = "enrichment",
                              sce_layer = spe_pb) 

top_DEGs$test <- gsub("_", "~", top_DEGs$test)

write.csv(top_DEGs, file = here(data_dir, "enrichment_modeling_SpD_top100.csv"), row.names = FALSE)

# slurmjobs::job_single('20_model_pseudobulk_anno', create_shell = TRUE, memory = '100G', command = "Rscript 20_model_pseudobulk_anno.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
        
        
