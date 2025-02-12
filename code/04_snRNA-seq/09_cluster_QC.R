## Louise Huuki-Myers, January 2025
## Run spatialLIBD::registration_wrapper to get cell type modeling and pseudobulk data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("tidyverse")
library("HDF5Array")
library("getopt")

## Set up dirs
data_dir <- here("processed-data", "04_snRNA-seq", "09_cluster_QC")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "09_cluster_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

## load cluster outputs
load(here("processed-data", "04_snRNA-seq", "11_cluster_sn", "walktrap_snn_k10_clusters.Rdata"), verbose = TRUE) ## temp
# load(here("processed-data", "04_snRNA-seq", "07_cluster_sn_prelim", "walktrap_snn_k%02d_clusters_prelim.Rdata"), verbose = TRUE)
sce$snn_k10 <- sprintf("k10c%02d", clusters)

#### scType annotations ####
load(here("processed-data", "04_snRNA-seq", "08_sctype_prelim", "sctype_prelim.Rdata"), verbose = TRUE)

# anno_table <- sctype[match(sce$snn_k10, sctype$cluster), c("cell_type_class", "cell_type_broad", "cell_type_fine")]

pd <- as.data.frame(colData(sce)) |>
    left_join(sctype |> select(snn_k10 = cluster, cell_type_class, cell_type_broad, cell_type_fine))


#### summary by cluster ####
message(Sys.time(), " - Summarize cluster")

cluster_info <- pd |>
    group_by(snn_k10, cell_type_class, cell_type_broad, cell_type_fine) |>
    summarize(n = n(),
              prop = n/ncol(sce),
              # cProp_APOE_carrier = sum(APOE_carrier == "E4+")/n,
              # cProp_Anc_AA = sum(Ancestry == "AA")/n,
              # cProp_Sex_F = sum(Sex == "F")/n,
              median_sum = median(sum),
              median_detected = median(detected),
              median_Mito_percent = median(subsets_Mito_percent),
              median_scDblFinder.score = median(scDblFinder.score)
              )

# check cutoffs 

cell_class_cutoffs <- pd |> 
    select(cell_type_class, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!cell_type_class, names_to = "metric") |>
    group_by(cell_type_class, metric) |>
    summarize(median = median(value),
              mad = mad(value)) |>
    mutate(cutoff = case_when(metric == "scDblFinder.score" ~ 0.5,
                              metric == "subsets_Mito_percent" ~ median + 3*mad,
                              TRUE ~ median - 1*mad))

cluster_metrics_long <- pd |> 
    select(cell_type_class, cell_type_broad, cell_type_fine, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(cell_type_fine,cell_type_broad, cell_type_class), names_to = "metric") |>
    group_by(cell_type_fine, cell_type_class, cell_type_broad, metric) |>
    summarize(median = median(value)) |>
    left_join(cell_class_cutoffs |> select(cell_type_class, metric, cutoff)) |>
    mutate(pass_metricQC = case_when(metric %in% c("scDblFinder.score", "subsets_Mito_percent") ~ median < cutoff,
                              TRUE ~ median > cutoff)) |>
    group_by(cell_type_fine, cell_type_class, cell_type_broad) |>
    mutate(passALL_metricQC = all(pass_metricQC))

cluster_metrics_long |> ungroup() |> count(pass_metricQC, passALL_metricQC)

## add pass fail to pd & cluster_info
cluster_metrics_pass <- cluster_metrics_long |> ungroup() |> select(cell_type_fine, passALL_metricQC)  |> unique()

pd <- pd |> left_join(cluster_metrics_pass)
cluster_info <- cluster_info |> left_join(cluster_metrics_pass)

cluster_info |>
    group_by(passALL_metricQC) |>
    summarize(n = sum(n),
              prop = n/nrow(pd))

cluster_info |>
    group_by(passALL_metricQC, cell_type_broad) |>
    summarize(n = sum(n))

write.csv(cluster_info, file = here(data_dir, "sn_clusterQC_info.csv"), row.names = FALSE)

#### Plot Quality Metrics ####

## median boxplots 

qc_median_boxplot <- cluster_metrics_long |>
    ggplot(aes(x = cell_type_broad, y = median)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(fill = cell_type_fine, shape = passALL_metricQC)) +
    geom_hline(aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    facet_grid(metric~cell_type_class, scales = "free") +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24)) +
    theme_bw()
    
ggsave(qc_median_boxplot, filename = here(plot_dir, "sn_clusterQC-QCmedian_boxplot.png"))    


## violin plots
qc_violins <- map(c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"), function(met){
    
    qc_violin_plot <- pd |>
        ggplot(aes(x = cell_type_fine, y = !!sym(met), fill = cell_type_fine, colour = passALL_metricQC)) +
        geom_violin(scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
        scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
        scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
        theme_bw() +
        geom_hline(data = cell_class_cutoffs |> filter(metric == met), aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
        facet_wrap(~cell_type_class, ncol = 2, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(qc_violin_plot, filename = here(plot_dir, sprintf("sn_clusterQC-QCmetricViolin_%s.png", met)), width = 12, height = 5)
    
    return(qc_violin_plot)
    
})


map(c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"), function(m){
   cell_class_cutoffs |> filter(metric == m) 
        
})

#### plot n cells and proportions ####

## colors 
cell_type_colors <- DeconvoBuddies::create_cell_colors(cell_types = levels(sctype$cell_type_fine), split = "\\.")

## n cell barplot
n_cell_barplot <- cluster_info |>
    ggplot(aes(x = cell_type_fine, y = n, fill = cell_type_fine, color = passALL_metricQC)) +
    geom_col() +
    geom_text(aes(label = n), size = 2, vjust = -0.5) +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(n_cell_barplot, filename = here(plot_dir, "sn_clusterQC-n_cells_barplot.png"), width = 10)

## proportion barplot


#### plot reduced dims ####
colData(sce) <- DataFrame(pd)

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

walk2(c("cell_type_broad", "cell_type_fine"), cell_type_colors, function(ct, colors){
    
    walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn_clusterQC",
                                                  var_type = "cat", 
                                                  dimred = .x, 
                                                  my_var = ct, 
                                                  color_pal = colors,
                                                  sufix = "prelim_k10"))
})



walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn_clusterQC",
                                              var_type = "cat", 
                                              dimred = .x, 
                                              my_var = "passALL_metricQC", 
                                              color_pal = c(`TRUE` = "black", `FALSE` = "red"),
                                              sufix = "prelim_k10"))

## filtered 
sce <- sce[,sce$passALL_metricQC]

walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn_clusterQCpass",
                                              var_type = "cat", 
                                              dimred = .x, 
                                              my_var = "cell_type_broad", 
                                              color_pal = cell_type_colors$broad,
                                              sufix = "prelim_k10"))


# slurmjobs::job_single('09_cluster_QC',
#                       create_shell = TRUE, memory = '10G', 
#                       command = "Rscript 09_cluster_QC.R --cluster 'ct_broad_k20'")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

