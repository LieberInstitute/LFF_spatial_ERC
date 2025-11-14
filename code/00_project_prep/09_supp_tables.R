## Louise Huuki-Myers, Oct 2025
## Compile format and export supp tables

library("tidyverse")
library("writexl")
library("here")
library("sessioninfo")
library(spatialLIBD)

data_dir <- here("processed-data", "00_project_prep", "09_supp_tables")
if(!dir.exists(data_dir)) dir.create(data_dir)

supp_tables <- list()

#### Table S1 Donor demographics. ####
supp_tables$TableS1 <- read.csv(here("processed-data", "00_project_prep", "05_pathology", "sample_taupathy.csv"))


#### Table S2 SRT sample information ####

supp_tables$TableS2 <- read.csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "ERC_Visium_summary_sample.csv"), row.names = 1)


#### Table S3 SpD Marker Genes ####

## SpD markers
load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_SpD.Rdata"), verbose = TRUE)
# marker_stats

erc_SpD_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))

get_sig_genes



## snRNA-seq sub-cluster metrics
"processed-data/04_snRNA-seq/28_subcluster_update_sce/ERC_sn_subcluster_info.csv"

## visium info
read.csv(here("processed-data", "02_build_spe", "sample_info.csv"))


## Stable_sn_crumblr snRNA-seq crumblr results 

## stable_DEG_Visium_GO Visium GO results 

## stable_DEG_sn_GO snRNA-seq GO results - broad

## stable_DEG_sn_fine_GO snRNA-seq GO results - fine

#### Write table ####
map(supp_tables, dim)

write_xlsx(supp_tables, path =  here(data_dir, "SuppTables.xlsx"))
