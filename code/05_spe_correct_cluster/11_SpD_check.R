## September 2024, Louise Huuki-Myers
## Check quality metrics and marker genes expression of SpatialDomains

library("spatialLIBD")
library("scater")
library("tidyverse")
library("DeconvoBuddies")
library("here")
library("sessioninfo")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### specify and load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))



## create color pallet
k <- 9

SpD_colors <- Polychrome::palette36.colors(k)
names(SpD_colors) <- sort(unique(spe$BayesSpace_PCA_Harmony_k09))

# Sp09D01   Sp09D02   Sp09D03   Sp09D04   Sp09D05   Sp09D06   Sp09D07   Sp09D08   Sp09D09 
# "#5A5156" "#E4E1E3" "#F6222E" "#FE00FA" "#16FF32" "#3283FE" "#FEAF16" "#B00068" "#1CFFCE" 

#### Check QC metrics by SpD ####

message(Sys.time(), "- Scran Outlier Violin Plots")

pdf(here(plot_dir, "spe_erc_QC_BayesSpace_PCA_Harmony_k09.pdf"), width = 21)

plotColData(spe, x = "BayesSpace_PCA_Harmony_k09", y = "sum_umi", colour_by = "scran_low_lib_size") +
    ggtitle("sum_umi") +
    scale_fill_manual(values = SpD_colors)

plotColData(spe, x = "BayesSpace_PCA_Harmony_k09", y = "sum_gene", colour_by = "scran_low_n_features") +
    ggtitle("sum_gene") +
    scale_fill_manual(values = SpD_colors)

plotColData(spe, x = "BayesSpace_PCA_Harmony_k09", y = "expr_chrM_ratio", colour_by = "scran_high_Mito_percent") +
    ggtitle("expr_chrM_ratio") +
    scale_fill_manual(values = SpD_colors) 

dev.off()

pd <- as.data.frame(colData(spe))

pd_qc_long <- pd |>
    select(sample_id, key, BayesSpace_PCA_Harmony_k09, sum_umi, sum_gene, expr_chrM_ratio) |>
    pivot_longer(!c("sample_id", "key", "BayesSpace_PCA_Harmony_k09"))

SpD_qc <- pd_qc_long |>
    ggplot(aes(x = BayesSpace_PCA_Harmony_k09, 
                             y = value, 
                             fill = BayesSpace_PCA_Harmony_k09)) +
    geom_boxplot() +
    facet_wrap(~name, ncol = 1, scales= "free_y") +
    scale_fill_manual(values = SpD_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(SpD_qc, filename = here(plot_dir, "SpD09_QC.png"))

#### Check marker genes by SpD ####

rownames(spe) <- rowData(spe)$gene_name

spe_pb_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_PCA_Harmony","spe_pseudobulk-BayesSpace_PCA_Harmony_k09.rds"))
modeling_results_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_PCA_Harmony","modeling_results-BayesSpace_PCA_Harmony_k09.rds"))

modeling_results_k09$enrichment |> 
    select(gene, ends_with("Sp09D07")) |>
    filter(logFC_Sp09D07 > 0) |>
    arrange(fdr_Sp09D07) |>
    head(10)

#                     gene t_stat_Sp09D07 p_value_Sp09D07  fdr_Sp09D07 logFC_Sp09D07
# ENSG00000197971      MBP       5.439940    1.161905e-07 1.391170e-07     1.5887006 * Oligo
# ENSG00000168314     MOBP       5.361880    1.723183e-07 2.049026e-07     1.6362505 * Oligo
# ENSG00000064787    BCAS1       5.139315    5.176398e-07 6.039571e-07     1.3500216
# ENSG00000136541     ERMN       4.772775    2.928887e-06 3.334874e-06     1.4077814
# ENSG00000160781    PAQR6       4.699450    4.093380e-06 4.637675e-06     0.9957585
# ENSG00000078804 TP53INP2       4.491110    1.036391e-05 1.159142e-05     0.9421917
# ENSG00000108387  SEPTIN4       4.427812    1.365368e-05 1.523207e-05     0.9601562
# ENSG00000197430   OPALIN       4.300728    2.352742e-05 2.610959e-05     1.2275028 * Oligodendrocytic Myelin Paranodal And Inner Loop Protein
# ENSG00000123560     PLP1       4.090538    5.629261e-05 6.195332e-05     1.2083869 * Oligo
# ENSG00000064393    HIPK2       4.047131    6.711421e-05 7.364323e-05     0.8292675

spe_pb_k09$spatialLIBD <- spe_pb_k09$BayesSpace_PCA_Harmony_k09

sig_genes <- sig_genes_extract_all(
    modeling_results = modeling_results_k09,
    sce_layer = spe_pb_k09,
    n = 10
) |>
    as.data.frame()

enrich_top <- sig_genes |>
    filter(model_type  == "enrichment", top <= 5)

enrich_top_list <- map(rafalib::splitit(enrich_top5$test), ~enrich_top5$gene[.x])
    
plot_marker_express_List(
    sce = spe,
    gene_list = enrich_top_list,
    pdf_fn = here(plot_dir, "SpD09_enrichment_top5.pdf"),
    cellType_col = "BayesSpace_PCA_Harmony_k09",
    gene_name_col = "gene_name",
    color_pal = SpD_colors,
    plot_points = FALSE
)

## from litature
modeling_results <- fetch_data(type = "modeling_results")

layer_marker_genes <- read.csv(here(data_dir, "layer_marker_genes.csv"))

plot_marker_express_List(
    sce = spe,
    gene_list = list(`Maynard_NatNeuro` = layer_marker_genes$gene),
    pdf_fn = here(plot_dir, "SpD09_layer_markers.pdf"),
    cellType_col = "BayesSpace_PCA_Harmony_k09",
    gene_name_col = "gene_name",
    color_pal = SpD_colors,
    plot_points = FALSE
)

sample_order <- sort(unique(spe$BrNum))
walk(c("AQP4", "HPCAL1", "KRT17", "MOBP", "PCP4", "SNAP25"),
          ~vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, sprintf("spe_erc-%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order)

     )




# slurmjobs::job_single('05_BayesSpace', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
