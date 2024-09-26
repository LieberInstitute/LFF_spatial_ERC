## September 2024, Louise Huuki-Myers
## Check quality metrics and marker genes expression of SpatialDomains

library("spatialLIBD")
library("scater")
library("tidyverse")
library("DeconvoBuddies")
library("bluster")
library("here")
library("sessioninfo")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### source plotting function ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

#### specify and load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

reducedDimNames(spe)
# [1] "10x_pca"      "10x_tsne"     "10x_umap"     "PCA_p1"       "PCA_p2"       "TSNE"         "UMAP"        
# [8] "HARMONY"      "UMAP.HARMONY" "TSNE.HARMONY"

## create color pallet

SpD_colors <- map(c(2:9), function(k){
    
    colors <- Polychrome::palette36.colors(k)
    if(k < 3) colors <- colors[1:2]
    names(colors) <- sort(unique(spe[[sprintf("BayesSpace_PCA_Harmony_k%02d", k)]]))
    return(colors)
})

names(SpD_colors) <- sprintf("k%02d", 2:9)



#### Plot clusters in reduced dims ####

walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){
my_plot_reduced_dim(spe, 
                    prefix = "spe", 
                    var_type = "cat", 
                    dimred = dim, 
                    my_var = "BayesSpace_PCA_Harmony_k02", 
                    color_pal = SpD_colors[1:2])
    })

walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){

    walk(c(2:9), function(k){
        
        SpD_colors_k <- SpD_colors[[sprintf("k%02d", k)]]
        
        my_plot_reduced_dim(spe, 
                            prefix = "spe",
                            var_type = "cat", 
                            dimred = dim, 
                            my_var = sprintf("BayesSpace_PCA_Harmony_k%02d", k), 
                            color_pal = SpD_colors_k)
    
        })
})
    


#### Check QC metrics by SpD ####

message(Sys.time(), "- Scran Outlier Violin Plots")

map(c(2,9), function(k){
    
    SpD_colors_k <- SpD_colors[[sprintf("k%02d", k)]]
    my_var = sprintf("BayesSpace_PCA_Harmony_k%02d", k)
    
    pdf(here(plot_dir, sprintf("spe_erc_QC_%s.pdf", my_var)), width = 21)
    
    plotColData(spe, x = my_var, y = "sum_umi", colour_by = "scran_low_lib_size") +
        ggtitle("sum_umi") +
        scale_fill_manual(values = SpD_colors_k)
    
    plotColData(spe, x = my_var, y = "sum_gene", colour_by = "scran_low_n_features") +
        ggtitle("sum_gene") +
        scale_fill_manual(values = SpD_colors_k)
    
    plotColData(spe, x = my_var, y = "expr_chrM_ratio", colour_by = "scran_high_Mito_percent") +
        ggtitle("expr_chrM_ratio") +
        scale_fill_manual(values = SpD_colors_k) 
    
    dev.off()
})



pd <- as.data.frame(colData(spe))

pd_qc_long <- map(c(k02 = 2, k09 = 9), function(k){
    
    SpD_colors_k <- SpD_colors[[sprintf("k%02d", k)]]
    my_var = sprintf("BayesSpace_PCA_Harmony_k%02d", k)
    
    pd_qc_long <- pd |>
        select(sample_id, key, sum_umi, sum_gene, expr_chrM_ratio, all_of(!!my_var)) |>
        pivot_longer(!c("sample_id", "key", my_var))
    
    my_var <- sym(my_var)
    
    SpD_qc <- pd_qc_long |>
        ggplot(aes(x = !!my_var, 
                   y = value, 
                   fill = !!my_var)
               ) +
        geom_boxplot() +
        facet_wrap(~name, ncol = 1, scales= "free_y") +
        scale_fill_manual(values = SpD_colors_k) +
        theme_bw() +
        theme(legend.position = "None")
    
    ggsave(SpD_qc, filename = here(plot_dir, sprintf("SpD%02d_QC.png", k)))
    
    return(pd_qc_long)
})



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

#### Jaccard indicies #####

table(spe$scran_discard, spe$BayesSpace_PCA_Harmony_k09)

## 70% of Sp09D07 are low library size
table(spe$scran_low_lib_size, spe$BayesSpace_PCA_Harmony_k09)
#       Sp09D01 Sp09D02 Sp09D03 Sp09D04 Sp09D05 Sp09D06 Sp09D07 Sp09D08 Sp09D09
# TRUE        2     329      18      72      15     114    3693       5    1669
# FALSE    8310    5666    8942   54805   10210   10086    1544    9697    7025

table(spe$BayesSpace_PCA_Harmony_k02, spe$BayesSpace_PCA_Harmony_k09)
#         Sp09D01 Sp09D02 Sp09D03 Sp09D04 Sp09D05 Sp09D06 Sp09D07 Sp09D08 Sp09D09
# Sp02D01    8312    5991    8960   54877   10225   10200    1720    9702    8648
# Sp02D02       0       4       0       0       0       0    3517       0      46

jacc.mat <- linkClustersMatrix(spe$BayesSpace_PCA_Harmony_k02, spe$BayesSpace_PCA_Harmony_k09)


# slurmjobs::job_single('05_BayesSpace', create_shell = TRUE, memory = '25G', command = "Rscript 05_BayesSpace.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
