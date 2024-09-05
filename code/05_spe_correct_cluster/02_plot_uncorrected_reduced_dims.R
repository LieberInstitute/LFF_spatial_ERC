# Aug, 2024 - Louise Huuki-Myers
# Plot uncorrected reduced dim plots for Visium data

library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scater")
library("HDF5Array")
library("viridis")
library("purrr")

plot_dir <- here("plots", "05_spe_correct_cluster", "02_plot_reduced_dims")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "02_build_spe", "01_preprocess_Harmony")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# spe <- readRDS(here("processed-data", "02_build_spe", "spe.rds"))

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
qc_colors <- c(fold = "#4BA402", out_edge = "#097FE0", vessel ="#BB1EF4", None = "#CCCCCC40")
scran_colors <- c(`TRUE` = "red", `FALSE` = "#CCCCCC40")

## define custom plotting function
my_plots <- my_plot_reduced_dim(spe = spe)

## categorical
walk(c("sample_id", "round", "APOE"), ~my_plot_reduced_dim(spe, dimred = "UMAP", my_var = .x, sufix = "uncorrected"))
walk(c("sample_id", "round", "APOE"), ~my_plot_reduced_dim(spe, dimred = "TSNE", my_var = .x, sufix = "uncorrected"))

## categorical + color scheme
walk2(c("qc_anno", "scran_discard"), list(qc_colors, scran_colors), ~my_plot_reduced_dim(spe, dimred = "UMAP", my_var = .x, sufix = "uncorrected", color_pal = .y))

## continuous
walk(c("sum_umi", "MBP"), ~my_plot_reduced_dim(spe, dimred = "UMAP", cat_var = FALSE, my_var = .x, sufix = "uncorrected"))

## expression

## plot layer marker genes
purrr::walk(c("MBP", "PCP4", "RELN", "SNAP25"), function(gene){
    umap_gene <- ggcells(spe, mapping = aes(x = UMAP.1, y = UMAP.2, color = !!sym(gene)))+
        geom_point(size = 0.2, alpha = 0.3) +
        coord_equal() +
        theme_bw() +
        scale_color_viridis(name = "logcounts") +
        labs(title = gene, x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
        theme(plot.title = element_text(face = "italic"))
    
    ggsave(umap_gene , filename = here(plot_dir, paste0("UMAP_uncorrected_expression_",gene,".png")))
})


# slurmjobs::job_single('02_plot_reduced_dims', create_shell = TRUE, memory = '25G', command = "Rscript 02_plot_reduced_dims.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

