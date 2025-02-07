## Louise Huuki-Myers, January 2025
## Run spatialLIBD::registration_wrapper to get cell type modeling and pseudobulk data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("tidyverse")
library("HDF5Array")
library("getopt")

# Import command-line parameters
# scec <- matrix(
#     c(  "cluster", "c", "1", "character", "Name of cluster"),
#     ncol = 5, byrow = TRUE
# )
# opt <- getopt(scec)

cluster <- opt$cluster

# cluster = "ct_fine_k20"

data_dir <- here("processed-data", "04_snRNA-seq", "13_cluster_check")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "13_cluster_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

pd <- as.data.frame(colData(sce))

#### summary by cluster ####
message(Sys.time(), " - Summarize cluster")
cluster_info <- pd |>
    group_by(!!sym(cluster)) |>
    summarize(n = n(),
              prop = n/ncol(sce),
              cProp_APOE_carrier = sum(APOE_carrier == "E4+")/n,
              cProp_Anc_AA = sum(Ancestry == "AA")/n,
              cProp_Sex_F = sum(Sex == "F")/n,
              median_sum = median(sum),
              median_detected = median(detected),
              median_Mito_percent = median(subsets_Mito_percent),
              median_scDbl_score = median(scDblFinder.score)
              )

#### plot n cells and proportions ####

cluster_colors <- metadata(sce)$cell_type_colors$fine

n_cell_barplot <- cluster_info |>
    ggplot(aes(x = !!sym(cluster), y = n, fill = !!sym(cluster))) +
    geom_col() +
    scale_fill_manual(values = cluster_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 
    
ggsave(n_cell_barplot, filename = here(plot_dir, sprintf("sn_%s-n_cells_barplot.png", cluster)))

#### Plot Quality Metrics ####

walk(c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"), function(metric){
    qc_violin_plot <- pd |>
        ggplot(aes(x = !!sym(cluster), y = !!sym(metric), fill = !!sym(cluster))) +
        geom_violin(scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
        scale_fill_manual(values = cluster_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "None") 
    
    ggsave(qc_violin_plot, filename = here(plot_dir, sprintf("sn_%s-QCviolin_%s.png", cluster, metric)), width = 11)
    
})



# slurmjobs::job_single('11_sn_model_pseudobulk', create_shell = TRUE, memory = '100G', command = "Rscript 11_sn_model_pseudobulk.R -cluster 'ct_broad_k20'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

