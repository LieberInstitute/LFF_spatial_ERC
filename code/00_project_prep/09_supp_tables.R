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

## visium info
read.csv(here("processed-data", "02_build_spe", "sample_info.csv"))

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
# read.csv(here("processed-data", "04_snRNA-seq", "33_sn_subcluster_summary", "ERC_sn_subcluster_summary_cell_type.csv"), row.names = 1) 

#### Table S6 Cell type subcluster marker genes ####

## cell_type_anno markers
load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio", "marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)
# marker_stats_MeanRatio

erc_sn_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_anno.rds"))
sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_pseudobulk-cell_type_anno.rds"))

sig_genes <- sig_genes_extract(
    modeling_results = erc_sn_modeling["enrichment"],
    sce_layer = sce_pb,
    n = nrow(sce_pb),
    gene_name = "gene_name",
)

# marker_stats_MeanRatio |> filter(gene_name == "MBP")
# sig_genes |> filter(ensembl == "ENSG00000197971")
# marker_stats |> filter(gene_ensembl == "ENSG00000197971")

head(sig_genes)
dim(sig_genes)
# [1] 941640     10 TOO BIG keep in sperate file

cell_type_markers <- sig_genes |> 
    select(cellType.target = test, ensembl, gene_name = gene, enrichment_stat = stat, enrichment_pval = pval, enrichment_fdr= fdr, enrichment_logFC = logFC) |>
    left_join(marker_stats_MeanRatio |>
                  select(cellType.target = cellType.target, ensembl = gene_ensembl, gene_name, MeanRatio, MeanRatio.rank, mean.target, cellType.2nd = cellType.2nd, mean.2nd) |>
                  mutate(cellType.target = as.character(cellType.target))) |>
    as_tibble()


cell_type_markers |> filter(gene_name == "MBP")

# write_csv(cell_type_markers, "/Users/louise.huuki/Library/CloudStorage/OneDrive-LieberInstituteforBrainDevelopment/LFF_ERC_paper/SuppTables/SuppTable6_cell_type_marker.csv")
# write_xlsx(cell_type_markers, "/Users/louise.huuki/Library/CloudStorage/OneDrive-LieberInstituteforBrainDevelopment/LFF_ERC_paper/SuppTables/SuppTable6_cell_type_marker.xlsx")

supp_tables$TableS6 <- cell_type_markers

#### Table S7 Spot deconvolution results ####
##TODO
list.files(here("processed-data", "15_spot_deconvolution", "03_results_RCTD", "cell_type_broad"))

#### Table S8 snRNA-seq cell type differential proportions. ####

supp_tables$TableS8 <- map_dfr(c("APOE_carrier", "Sex", "Age"), 
    ~read.csv(here("processed-data", "04_snRNA-seq", "22_crumblr_sn", paste0("sn_diff_prop_tree_test_",.x,".csv"))) |>
        mutate(test = .x)
    )

#### Table S9, 11, 13 differential expression results. ####
datatypes <- c("Visium", "sn_broad", "sn_fine")
    
DEG_carrier <- map(datatypes, ~readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", .x, sprintf("DGE_results_carrier_%s.Rds", .x))) |> 
                       dplyr::select(1:20) |>
                       dplyr::rename_with(~gsub("vlmf", "carrier", .))
                   )

DEG_ancestry <- map(datatypes, ~readRDS(here("processed-data", "13_compile_DGE", "05_compile_DGE_ancestry", .x, sprintf("DGE_results_ancestry_%s.Rds", .x))) |> 
                        dplyr::select(1:20) |>
                        dplyr::rename_with(~gsub("vlmf", "ancestry", .))
                    )

DEG_Sex <- map(datatypes, ~readRDS(here("processed-data", "13_compile_DGE", "09_compile_DGE_Sex", .x, sprintf("DGE_results_Sex_%s.Rds", .x))) |> 
                   dplyr::select(1:20) |>
                   dplyr::rename_with(~gsub("vlmf", "sex", .))
               )

combined_data <- pmap(list(carrier = DEG_carrier, anc = DEG_ancestry, sex = DEG_Sex), function(carrier, anc, sex){
    carrier |>
        left_join(anc) |>
        left_join(sex)
})

map(combined_data, ncol)

supp_tables$TableS9 <- combined_data[[1]]
supp_tables$TableS11 <- combined_data[[2]]
supp_tables$TableS13 <- combined_data[[3]]

rm(combined_data)

#### Table S10, 12, 14 GO over-representation analysis. ####

GO_carrier <- 

#### Write table ####

supp_tables <- supp_tables[sort(names(supp_tables))]

map(supp_tables, dim)

write_xlsx(supp_tables, path =  "/Users/louise.huuki/Library/CloudStorage/OneDrive-LieberInstituteforBrainDevelopment/LFF_ERC_paper/SuppTables/SuppTables.xlsx")
