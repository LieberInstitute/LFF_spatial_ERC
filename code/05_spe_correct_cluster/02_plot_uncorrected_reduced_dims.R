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
plot_dir <- here("plots", "05_spe_correct_cluster", "02_plot_uncorrected_reduced_dims")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "02_build_spe", "01_preprocess_Harmony")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))
spe
# class: SpatialExperiment 
# dim: 30494 122202 
# metadata(0):
#     assays(2): counts logcounts
# rownames(30494): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(122202): AAACAACGAATAGTTC-1_Br5212 AAACAAGTATCTCCCA-1_Br5212 ... TTGTTTGTATTACACG-1_Br6263
# TTGTTTGTGTAAATTC-1_Br6263
# colData names(48): sample_id in_tissue ... scran_quick_cluster sizeFactor
# reducedDimNames(7): 10x_pca 10x_tsne ... TSNE UMAP
# mainExpName: NULL
# altExpNames(0):
#     spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

reducedDimNames(spe)
# [1] "10x_pca"  "10x_tsne" "10x_umap" "PCA_p1"   "PCA_p2"   "TSNE"     "UMAP"  

## swap rownames to gene names for plotting markers
rownames(spe) <- rowData(spe)$gene_name

#### Establish color schemes ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

qc_colors <- c(fold = "#4BA402", out_edge = "#097FE0", vessel ="#BB1EF4", None = "#CCCCCC40")
scran_colors <- c(`TRUE` = "red", `FALSE` = "#CCCCCC40")


walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){
    
    sufix = "PCA_p1"
    
    ## categorical
    walk(c("sample_id", "round", "Visium_slide"), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = sufix))

    ## categorical facet
    walk(c("sample_id", "round", "Visium_slide"), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = sufix, facet = TRUE))

    ## categorical + color scheme
    walk2(c("qc_anno", "scran_discard", "APOE", "Ancestry", "Sex"),
          list(qc_colors, scran_colors, APOE_genotype_colors, ancestry_colors, sex_colors),
          ~my_plot_reduced_dim(spe,prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = sufix, color_pal = .y))

    walk2(c("APOE", "Ancestry", "Sex"),
          list(APOE_genotype_colors, ancestry_colors, sex_colors),
          ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = sufix, facet = TRUE, color_pal = .y))

    ## continuous
    walk(c("Age", "Rin"),
         ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "con", dimred = dim, my_var = .x, sufix = sufix))


    # ## continuous QC metrics
    walk(c("sum_umi", "sum_gene", "expr_chrM_ratio"),
         ~my_plot_reduced_dim(spe, prefix = "spe_QC", var_type = "con", dimred = dim, my_var = .x, sufix = sufix))

    # ## plot layer marker genes
    purrr::walk(c("MBP", "PCP4", "RELN", "SNAP25"), function(gene){
        my_plot_reduced_dim(spe, prefix = "spe_gene", var_type = "express", dimred = dim, my_var = gene, sufix = sufix)
    })
})


annotations <- unique(map_dfr(list.files(here("processed-data", "05_spe_correct_cluster","reduced_dim_anno"), full.names = TRUE), read.csv))

annotations |> dplyr::count(sample_id, ManualAnnotation)
# sample_id ManualAnnotation    n
# 1    Br1691    UMAP_outliers 3109
# 2    Br5276   UMAP_high_mito 1975
# 3    Br5529   UMAP_high_mito 1692
# 4    Br5599    UMAP_outliers 4209
# 5    Br5634    UMAP_outliers 3442
# 6    Br6085    UMAP_outliers   69
# 7    Br6538    UMAP_outliers 1818

annotations2 <- annotations |> 
    mutate(ManualAnnotation_BrNum = paste(ManualAnnotation, sample_id)) |>
    as_tibble() |>
    column_to_rownames("spot_name") |>
    select(-sample_id)

annotations2 <- annotations2[colnames(spe), ]

colData(spe) <- cbind(colData(spe), annotations2)


walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){
    
    sufix = "PCA_p1"
    
    # ## categorical
    walk(c("ManualAnnotation", "ManualAnnotation_BrNum"), ~my_plot_reduced_dim(spe, prefix = "spe", var_type = "cat", dimred = dim, my_var = .x, sufix = sufix))

})


sample_ids <- sort(unique(spe$sample_id))
spatialLIBD::vis_grid_clus(
    spe = spe,
    clustervar = "ManualAnnotation",
    pdf_file = here(plot_dir, "UMAP_outlier_spot_plot-ALL.pdf"),
    sort_clust = FALSE,
    spatial = FALSE,
    point_size = 1,
    sample_order = sample_ids
)

# slurmjobs::job_single('02_plot_uncorrected_reduced_dims', create_shell = TRUE, memory = '25G', command = "Rscript 02_plot_uncorrected_reduced_dims.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

