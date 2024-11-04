# Nov, 2024 - Louise Huuki-Myers
# Add Use SVG + Marker genes to calc reduced dimesions run Harmony for BayesSpace input

library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scran")
library("scater")
library("harmony")
library("BiocParallel")
library("purrr")
library("scry")
library("HDF5Array")

plot_dir <- here("plots", "05_spe_correct_cluster", "17_SVG_PCA_Harmony")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "05_spe_correct_cluster", "17_SVG_PCA_PCA_Harmony")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

spe <- loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))
reducedDimNames(spe)

## clear Reduced Dims
spe$TSNE <- NULL
spe$UMAP <- NULL
spe$HARMONY <- NULL
spe$UMAP.HARMONY <- NULL
spe$TSNE.HARMONY <- NULL

#### load SVGs ####
load(here("processed-data", "05_spe_correct_cluster", "13_gather_nnSVG", "top.svg.plusMarkers.Rdata"), verbose = TRUE)

#### Reduced Dims ####
num_reduced_dims <- 50

# compute PCA
message(Sys.time(), " - Compute PCA")

## Run PCA with selected set of SVGs + marker genes
spe <- runPCA(spe, 
              subset_row = top.svg.plusMarkers,
              ncomponents = num_reduced_dims,
              name = "PCA_SVG")

## Plot PC1 vs. PC2
message(Sys.time(), " - PCA plots")
pca_sample <- plotReducedDim(spe, dimred = "PCA_p1",
                             color_by = "sample_id") +
    theme_bw()
ggsave(pca_sample, filename = here(plot_dir, "PCA_sample.png"), width = 10)

## plot PCA elbow
pca_names <- reducedDimNames(spe)[grep("PCA", reducedDimNames(spe))]

precent_var <- map_dfr(pca_names, ~data.frame(geneset = .x,
                                              PC = seq(ncol(reducedDim(spe, .x))),
                                              precentVar = attr(reducedDim(spe, .x), "percentVar")))

pca_elbow <- ggplot(precent_var, aes(PC, precentVar, color = geneset)) +
    geom_line() +
    theme_bw()

ggsave(pca_elbow, filename = here(plot_dir, "PCA_elbow.png"))

precent_var |>
    dplyr::group_by(geneset) |>
    dplyr::summarise(elbow = PCAtools::findElbowPoint(percent.var))

# geneset elbow
# <chr>   <int>
# 1 PCA_SVG     6
# 2 PCA_p1      6
# 3 PCA_p2      6

# percent.var <- attr(reducedDim(spe, "PCA_SVG"), "percentVar")
# chosen.elbow <- PCAtools::findElbowPoint(percent.var)
# message("PCA elbow: ", chosen.elbow)
# PCA elbow: 6

## TSNE & UMAP
message(Sys.time(), " - TSNE")
spe <-
    runTSNE(spe,
            dimred = "PCA_SVG",
            name = "TSNE"
    )

message(Sys.time(), " - UMAP")
spe <- runUMAP(spe, 
               dimred = "PCA_SVG")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
message(Sys.time(), " - Done Reduced Dims...")

#### Harmony Batch Correction ####

## Harmony uses the PCA slot - use SVG + marker PCA
reducedDim(spe, "PCA") <- reducedDim(spe, "PCA_SVG")

message(Sys.time(), " - Harmony")
spe <- RunHarmony(spe, "sample_id")

## ADD UMAP + TSNE
message(Sys.time(), " - Harmony UMAP")
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")

colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

message(Sys.time(), " - HARMONY TSNE")
spe <- runTSNE(spe,
               dimred = "HARMONY",
               name = "TSNE.HARMONY")

## Remove redundant PCA
reducedDim(spe, "PCA") <- NULL

#### Save data ####
message(Sys.time(), " - Saving HDF5 SPE")
saveHDF5SummarizedExperiment(
    spe,
    dir = here("processed-data", "spe_objects", "spe_PCA_SVG"),
    replace = TRUE
)

#### Plot Reduced Dims ####
## source plotting function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## swap rownames to gene names for plotting markers
rownames(spe) <- rowData(spe)$gene_name

#### Establish color schemes ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

qc_colors <- c(fold = "#4BA402", out_edge = "#097FE0", vessel ="#BB1EF4", None = "#CCCCCC40")
scran_colors <- c(`TRUE` = "red", `FALSE` = "#CCCCCC40")


walk(c("UMAP", "TSNE", "UMAP.HARMONY", "TSNE.HARMONY"), function(dim){
    
    sufix = "SVGm"
    
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
    
    
    ## continuous QC metrics
    walk(c("sum_umi", "sum_gene", "expr_chrM_ratio"),
         ~my_plot_reduced_dim(spe, prefix = "spe_QC", var_type = "con", dimred = dim, my_var = .x, sufix = sufix))
    
    ## plot layer marker genes
    purrr::walk(c("MBP", "PCP4", "RELN", "SNAP25"), function(gene){
        my_plot_reduced_dim(spe, prefix = "spe_gene", var_type = "express", dimred = dim, my_var = gene, sufix = sufix)
    })
})


# slurmjobs::job_single('17_SVG_PCA_Harmony', create_shell = TRUE, memory = '25G', command = "Rscript 17_SVG_PCA_Harmony.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
