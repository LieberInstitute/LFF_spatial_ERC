# Aug, 2024 - Louise Huuki-Myers
# Plot uncorrected reduced dim plots for Visium data

library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scater")
library("HDF5Array")
library("viridis")
library("purrr")

#### source plotting function ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

#### prep dirs ####
plot_dir <- here("plots", "05_spe_correct_cluster", "04_plot_harmony_reduced_dims")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_GLM_Harmony"))
spe

reducedDimNames(spe)
# [1] "10x_pca"          "10x_tsne"         "10x_umap"         "PCA_p1"           "PCA_p2"           "TSNE"            
# [7] "UMAP"             "GLMPCA_approx_1k" "GLMPCA_approx_2k" "GLMPCA_approx_5k" "HARMONY"          "UMAP.HARMONY"    
# [13] "TSNE.HARMONY"

## swap rownames to gene names for plotting markers
rownames(spe) <- rowData(spe)$gene_name

#### Establish color schemes ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

qc_colors <- c(fold = "#4BA402", out_edge = "#097FE0", vessel ="#BB1EF4", None = "#CCCCCC40")
scran_colors <- c(`TRUE` = "red", `FALSE` = "#CCCCCC40")


walk(c("UMAP.HARMONY", "TSNE.HARMONY"), function(dim){
    ## categorical
    walk(c("sample_id", "round", "Visium_slide"), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = "GLM"))
    
    ## categorical facet
    walk(c("sample_id", "round"), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = "GLM", facet = TRUE))
    walk2(c("APOE"), list(APOE_genotype_colors), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = "GLM", facet = TRUE, color_pal = .y))
    
    ## categorical + color scheme
    walk2(c("qc_anno", "scran_discard", "APOE", "Ancestry"),
          list(qc_colors, scran_colors, APOE_genotype_colors, ancestry_colors),
          ~my_plot_reduced_dim(spe,prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = "GLM", color_pal = .y))
    
    
    
    ## continuous QC metrics
    walk(c("sum_umi", "sum_gene", "expr_chrM_ratio"),
         ~my_plot_reduced_dim(spe, prefix = "spe_QC", var_type = "con", dimred = dim, my_var = .x, sufix = "GLM"))
    
    ## plot layer marker genes
    purrr::walk(c("MBP", "PCP4", "RELN", "SNAP25"), function(gene){
        my_plot_reduced_dim(spe, prefix = "spe_gene", var_type = "express", dimred = dim, my_var = gene, sufix = "GLM")
    })
})

#### UMAP annotations ####

head(reducedDims(spe)$UMAP.HARMONY)

outliers <- rownames(reducedDims(spe)$UMAP.HARMONY[reducedDims(spe)$UMAP.HARMONY[,"UMAP1"] < -5,])
length(outliers)
# [1] 1524

spe$GLM_outlier <- spe$key %in% outliers

table(spe$GLM_outlier, spe$sample_id)
table(spe$GLM_outlier, spe$scran_discard)
#         TRUE  FALSE
# FALSE   5284 115394
# TRUE    1345    179

sample_ids <- sort(unique(spe$sample_id))
spatialLIBD::vis_grid_clus(
    spe = spe,
    clustervar = "GLM_outlier",
    pdf_file = here(plot_dir, "GLM_outlier_spot_plot-ALL.pdf"),
    sort_clust = FALSE,
    colors = c(`TRUE` = "blue", `FALSE` = "#CCCCCC40"),
    spatial = FALSE,
    point_size = 1,
    sample_order = sample_ids
)

walk(c("UMAP.HARMONY", "TSNE.HARMONY"), function(dim){

    ## categorical + color scheme
    walk2(c("GLM_outlier"),
          list(c(`TRUE` = "blue", `FALSE` = "#CCCCCC40")),
          ~my_plot_reduced_dim(spe,prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = "GLM", color_pal = .y))
    })

# slurmjobs::job_single('04_plot_harmony_reduced_dims', create_shell = TRUE, memory = '25G', command = "Rscript 04_plot_harmony_reduced_dims.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


