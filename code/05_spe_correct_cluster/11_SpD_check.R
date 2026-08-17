## September 2024, Louise Huuki-Myers
## Check quality metrics and marker genes expression of Spatial Domains

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

SpD_colors <- map(c(2:28), function(k){
    
    colors <- Polychrome::palette36.colors(k)
    if(k < 3) colors <- colors[1:2]
    names(colors) <- sort(unique(spe[[sprintf("BayesSpace_SVGm_k%02d", k)]]))
    return(colors)
})

names(SpD_colors) <- sprintf("k%02d", 2:28)


table(spe$BayesSpace_SVGm_k02)

#### Plot clusters in reduced dims ####

# walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){
# my_plot_reduced_dim(spe, 
#                     prefix = "spe", 
#                     var_type = "cat", 
#                     dimred = dim, 
#                     my_var = "BayesSpace_SVGm_k02", 
#                     color_pal = SpD_colors$k02)
#     })

walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){

    walk(c(2:12), function(k){
        
        SpD_colors_k <- SpD_colors[[sprintf("k%02d", k)]]
        
        my_plot_reduced_dim(spe, 
                            prefix = "spe",
                            var_type = "cat", 
                            dimred = dim, 
                            my_var = sprintf("BayesSpace_SVGm_k%02d", k), 
                            color_pal = SpD_colors_k)
    
        })
})
    


#### Check QC metrics by SpD ####

message(Sys.time(), "- Scran Outlier Violin Plots")

map(c(2,12), function(k){
    
    SpD_colors_k <- SpD_colors[[sprintf("k%02d", k)]]
    my_var = sprintf("BayesSpace_SVGm_k%02d", k)
    
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

pd_qc_long <- map(c(k02 = 2, k09 = 9, k11 = 11), function(k){
    
    SpD_colors_k <- SpD_colors[[sprintf("k%02d", k)]]
    my_var = sprintf("BayesSpace_SVGm_k%02d", k)
    
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

spe_pb_k11 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k11.rds"))
spe_pb_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k09.rds"))
spe_pb_k02 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k02.rds"))

modeling_results_k11 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","modeling_results-BayesSpace_SVGm_k11.rds"))
modeling_results_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","modeling_results-BayesSpace_SVGm_k09.rds"))
modeling_results_k09 <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","modeling_results-BayesSpace_SVGm_k02.rds"))

modeling_results_k09$enrichment |> 
    select(gene, ends_with("Sp09D07")) |>
    filter(logFC_Sp09D07 > 0) |>
    arrange(fdr_Sp09D07) |>
    head(10)


spe_pb_k09$spatialLIBD <- spe_pb_k09$BayesSpace_SVGm_k09

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
    cellType_col = "BayesSpace_SVGm_k09",
    gene_name_col = "gene_name",
    color_pal = SpD_colors,
    plot_points = FALSE
)

sig_genes_k02 <- sig_genes_extract_all(
    modeling_results = modeling_results_k02,
    sce_layer = spe_pb_k02,
    n = 10
) |>
    as.data.frame()


modeling_results_k02$enrichment |>     
    dplyr::filter(logFC_Sp02D02 > 0) |>
    dplyr::arrange(fdr_Sp02D02) |>
    head()

modeling_results_k02$enrichment |>     
    dplyr::filter(logFC_Sp02D01 > 0) |>
    dplyr::arrange(fdr_Sp02D01) |>
    head()

# write_csv(modeling_results_k02$enrichment, file = here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "modeling_results-BayesSpace_SVGm_k02.csv"))

## from literature
modeling_results <- fetch_data(type = "modeling_results")

layer_marker_genes <- read.csv(here(data_dir, "layer_marker_genes.csv"))

plot_marker_express_List(
    sce = spe,
    gene_list = list(`Maynard_NatNeuro` = layer_marker_genes$gene),
    pdf_fn = here(plot_dir, "SpD09_layer_markers.pdf"),
    cellType_col = "BayesSpace_SVGm_k09",
    gene_name_col = "gene_name",
    color_pal = SpD_colors,
    plot_points = FALSE
)

sample_order <- sort(unique(spe$BrNum))

## SVGs

nnSVG_avg <- read.csv(here("processed-data", "05_spe_correct_cluster", "13_gather_nnSVG", "nnSVG_avg.csv"))

nnSVG_avg$gene_name

walk(top.svg,
     ~vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, sprintf("spe_erc-SVG_%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order))
     
#### position ####
spe_position <- read.table(here("processed-data", "02_build_spe", "spe_position.csv"), sep = "\t", header = TRUE)

# missing order
sample_order[!sample_order %in% spe_position$sample_id]
# [1] "Br1039" "Br1556" "Br1691" "Br1706" "Br2305" "Br3974" "Br5460" "Br5529" "Br5599" "Br5832" "Br5854" "Br5941"
# [13] "Br6085" "Br6098" "Br6476"

rownames(spe_position) <- spe_position$sample_id

spe$position <- spe_position[spe$sample_id, "position"]

spe$sample_id <- paste(spe$position, spe$sample_id)

table(spe$position)

dir.create(here(plot_dir, "position"))

sample_order <- sort(unique(spe$sample_id))
# walk(c("MBP", "AQP4", "HPCAL1", "KRT17", "MOBP", "PCP4", "SNAP25"),
walk(c("GULP1","COL25A1", "SATB1", "PEX5L"),
     ~vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, "position", sprintf("spe_erc_AMY-%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order)
     
)

walk(c("FN1", "NTS", "TLE1"),
     ~vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, "position", sprintf("spe_erc_HP-%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order)
     
)

#### spot plots ####
dir.create(here(plot_dir, "spot_plots"))

sample_ids <- sort(unique(spe$sample_id))

p_list = vis_grid_clus(
    spe = spe,
    clustervar = 'BayesSpace_SVGm_k02',
    pdf_file = here(plot_dir, "spot_plots", "spe_BayesSpace_SVGm_k02-ALL.pdf"),
    sort_clust = FALSE,
    spatial = FALSE,
    point_size = 1,
    sample_order = sample_ids
)


#### Jaccard indicies #####

table(spe$scran_discard, spe$BayesSpace_SVGm_k09)

## 70% of Sp09D07 are low library size
table(spe$scran_low_lib_size, spe$BayesSpace_SVGm_k09)
#       Sp09D01 Sp09D02 Sp09D03 Sp09D04 Sp09D05 Sp09D06 Sp09D07 Sp09D08 Sp09D09
# TRUE        2     329      18      72      15     114    3693       5    1669
# FALSE    8310    5666    8942   54805   10210   10086    1544    9697    7025

table(spe$BayesSpace_SVGm_k02, spe$BayesSpace_SVGm_k09)
#         Sp09D01 Sp09D02 Sp09D03 Sp09D04 Sp09D05 Sp09D06 Sp09D07 Sp09D08 Sp09D09
# Sp02D01    8312    5991    8960   54877   10225   10200    1720    9702    8648
# Sp02D02       0       4       0       0       0       0    3517       0      46

jacc.mat <- linkClustersMatrix(spe$BayesSpace_SVGm_k02, spe$BayesSpace_SVGm_k09)

# slurmjobs::job_single('011_SpD_check', create_shell = TRUE, memory = '25G', command = "Rscript 011_SpD_check.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
