## Dec 2024, Louise Huuki-Myers
## Extract and plot top enrichment model genes from selected cluster

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "05_spe_correct_cluster", "20_top_enrichment_genes")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "20_top_enrichment_genes")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load SPE data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))
# spe <- readRDS(here("code", "06_iSEE_app", "sce_ERC_iSEE.rds"))


#### Extract Top Layer Enrichment Genes ####


#### Top Enrichment Gene Heatmap ####


#### Spot plots ####



# slurmjobs::job_single('20_top_enrichment_genes', create_shell = TRUE, memory = '25G', command = "Rcript 20_top_enrichment_genes.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
