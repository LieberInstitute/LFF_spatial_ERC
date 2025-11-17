## Louise Huuki-Myers, Oct 2025
## Compile format and export supp tables

library("tidyverse")
library("writexl")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("SingleCellExperiment")

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
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "spe_pseudobulk-SpD.rds"))

sig_genes <- sig_genes_extract(
    modeling_results = erc_SpD_modeling["enrichment"],
    sce_layer = spe_pb,
    n = nrow(spe_pb),
    gene_name = "gene_id",
)

head(sig_genes)
dim(sig_genes)
# [1] 112284     10

# marker_stats |> filter(gene_name == "MBP")
# sig_genes |> filter(ensembl == "ENSG00000197971")
# marker_stats |> filter(gene_ensembl == "ENSG00000197971")

SpD_markers <- sig_genes |> 
    select(SpD.target = test, ensembl, enrichment_stat = stat, enrichment_pval = pval, enrichment_fdr= fdr, enrichment_logFC = logFC) |>
    mutate(SpD.target = gsub("_", "~", SpD.target)) |>
    left_join(marker_stats |> 
                  select(SpD.target = cellType.target, ensembl = gene_ensembl, gene_name, MeanRatio, MeanRatio.rank, mean.target, SpD.2nd = cellType.2nd,  mean.2nd) |> 
                  mutate(SpD.target = as.character(SpD.target))) |>
    as_tibble()

SpD_markers |> filter(!is.na(enrichment_stat), !is.na(MeanRatio))

SpD_markers |> filter(gene_name == "MBP")

supp_tables$TableS3 <- SpD_markers

#### Table S4 SpD Differential Proportions ####

supp_tables$TableS4 <- read.csv(here("processed-data", "05_spe_correct_cluster", "26_crumblr_SpD", "SpD_diff_prop_tree_test_APOE_carrier.csv"))

#### Table S5 snRNA-seq sample information ####
sn_sample_qc <-read.csv(here("processed-data", "04_snRNA-seq", "02_droplet_QC", "erc_sn_sample_QC_summary.csv"))
sn_sample_summary <- read.csv(here("processed-data", "04_snRNA-seq", "33_sn_subcluster_summary", "ERC_sn_subcluster_summary_sample.csv"), row.names = 1)  |>
    dplyr::rename(BrNum = sample_id)

supp_tables$TableS5 <- sn_sample_qc |> left_join(sn_sample_summary)

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
