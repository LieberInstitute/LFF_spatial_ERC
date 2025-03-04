## Louise Huuki-Myers, January 2025
## Run spatialLIBD::registration_wrapper to get cell type modeling and pseudobulk data

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("tidyverse")
library("HDF5Array")
library("ggrepel")

## Set up dirs
data_dir <- here("processed-data", "04_snRNA-seq", "09_cluster_QC")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "09_cluster_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

## load cluster outputs
load(here("processed-data", "04_snRNA-seq", "07_cluster_sn_prelim", "walktrap_snn_k10_clusters_prelim.Rdata"), verbose = TRUE)
# clusters
sce$snn_k10 <- sprintf("k10c%02d", clusters)

#### scType annotations ####
load(here("processed-data", "04_snRNA-seq", "08_sctype_prelim", "sctype_prelim.Rdata"), verbose = TRUE)
# sctype

table(sctype$cell_type_broad)
# Astro  Endo Micro Oligo   OPC Excit Inhib   Neu 
# 9     1     5     6     3    14     5     1 

pd <- as.data.frame(colData(sce)) |>
    left_join(sctype |> select(snn_k10 = cluster, cell_type_class, cell_type_broad, cell_type_fine))

pd |> filter(is.na(cell_type_fine))

all(sctype$cluster %in% unique(pd$snn_k10))
all(unique(pd$snn_k10) %in% sctype$cluster)

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
                              TRUE ~ median - 1*mad),
           cutoff_anno = case_when(metric %in% c("scDblFinder.score", "subsets_Mito_percent") ~  sprintf("< %g", round(cutoff, 2)),
                                   TRUE ~ sprintf("> %g", round(cutoff, 2))))

# cell_type_class metric                     median          mad   cutoff cutoff_anno
# <chr>           <chr>                       <dbl>        <dbl>    <dbl> <chr>      
# 1 glia            detected              2321         1103.       1218.    > 1217.95  
# 2 glia            scDblFinder.score        0.000996     0.000757    0.5   < 0.5      
# 3 glia            subsets_Mito_percent     0.0879       0.104       0.400 < 0.4      
# 4 glia            sum                   4957         3616.       1341.    > 1340.94  
# 5 neuron          detected              6711         2537.       4174.    > 4174.27  
# 6 neuron          scDblFinder.score        0.00197      0.00203     0.5   < 0.5      
# 7 neuron          subsets_Mito_percent     0.115        0.139       0.531 < 0.53     
# 8 neuron          sum                  28592        23000.       5592.    > 5592.43   

## save out cutoff data
write.csv(cell_class_cutoffs, file = here(data_dir,"ERC_sn_cell_class_cutoffs.csv"), row.names = FALSE)

cluster_metrics_long <- pd |> 
    select(snn_k10, cell_type_class, cell_type_broad, cell_type_fine, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(snn_k10, cell_type_fine,cell_type_broad, cell_type_class), names_to = "metric") |>
    group_by(snn_k10, cell_type_fine, cell_type_class, cell_type_broad, metric) |>
    summarize(median = median(value)) |>
    left_join(cell_class_cutoffs |> select(cell_type_class, metric, cutoff, cutoff_anno)) |>
    mutate(pass_metricQC = case_when(metric %in% c("scDblFinder.score", "subsets_Mito_percent") ~ median < cutoff,
                              TRUE ~ median > cutoff),
           metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))
           )|>
    group_by(snn_k10, cell_type_fine, cell_type_class, cell_type_broad) |>
    mutate(passALL_metricQC = all(pass_metricQC))

cluster_metrics_long |> ungroup() |> count(pass_metricQC, passALL_metricQC)

cluster_metrics_long |> ungroup() |> count(metric, cutoff_anno)

## add pass fail to pd & cluster_info
## TODO add list of failed metrics
cluster_metrics_pass <- cluster_metrics_long |> ungroup() |> select(cell_type_fine, passALL_metricQC)  |> unique()

table(cluster_metrics_pass$passALL_metricQC)

pd <- pd |> left_join(cluster_metrics_pass)
cluster_info <- cluster_info |> left_join(cluster_metrics_pass)

cluster_info |>
    group_by(passALL_metricQC) |>
    summarize(n = sum(n),
              prop = n/nrow(pd),
              n_clusters = length(unique(cell_type_fine)))

# passALL_metricQC      n  prop n_clusters
# <lgl>             <int> <dbl>      <int>
# 1 FALSE             14304 0.102         15
# 2 TRUE             125815 0.898         29

cluster_info |>
    group_by(passALL_metricQC, cell_type_broad) |>
    summarize(n = sum(n),
              n_clusters = length(unique(cell_type_fine)))

# passALL_metricQC cell_type_broad     n n_clusters
# <lgl>            <fct>           <int>      <int>
# 1 FALSE            Astro            1688          6
# 2 FALSE            Micro            4289          3
# 3 FALSE            Oligo            5127          3
# 4 FALSE            OPC               251          2
# 5 FALSE            Excit            2949          1
# 6 TRUE             Astro           22554          3
# 7 TRUE             Endo             1672          1
# 8 TRUE             Micro           13480          2
# 9 TRUE             Oligo           47203          3
# 10 TRUE             OPC             13680          1
# 11 TRUE             Excit           15975         13
# 12 TRUE             Inhib           10831          5
# 13 TRUE             Neu               420          1

#### Plot Quality Metrics ####
cell_type_colors <- DeconvoBuddies::create_cell_colors(levels(sctype$cell_type_fine), split = "\\.")

## quality scatter
median_sum_v_detected <- cluster_info |>
    ggplot(aes(median_sum, median_detected, fill = cell_type_fine, shape = passALL_metricQC)) +
    geom_point() +
    geom_text_repel(aes(label = cell_type_fine), size = 2) +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24)) +
    theme_bw()

ggsave(median_sum_v_detected, filename = here(plot_dir, "sn_clusterQC-scater_median_sum_v_detected.png"))


## median boxplots 

qc_median_boxplot <- cluster_metrics_long |>
    ggplot(aes(x = cell_type_broad, y = median)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(fill = cell_type_fine, shape = passALL_metricQC)) +
    geom_hline(aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    geom_text(aes(y = cutoff, label = cutoff_anno), x = 1.5, color = "blue", vjust = -.5) +
    facet_grid(metric~cell_type_class, scales = "free") +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24)) +
    theme_bw()
    
ggsave(qc_median_boxplot, filename = here(plot_dir, "sn_clusterQC-QCmedian_boxplot.png"))    


## violin plots
qc_violins <- map(c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"), function(met){
    
    qc_violin_plot <- pd |>
        ggplot() +
        geom_violin(aes(x = cell_type_fine, y = !!sym(met), fill = cell_type_fine, colour = passALL_metricQC),
                    scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
        scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
        scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
        theme_bw() +
        geom_hline(data = cell_class_cutoffs |> filter(metric == met), aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
        geom_label(data = cell_class_cutoffs |> filter(metric == met), aes(y = cutoff, label = cutoff_anno), x = 1.5, color = "blue", vjust = -.5) +
        facet_wrap(~cell_type_class, ncol = 2, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(qc_violin_plot, filename = here(plot_dir, sprintf("sn_clusterQC-QCmetricViolin_%s.png", met)), width = 12, height = 5)
    
    return(qc_violin_plot)
    
})

qc_violin_plot_all <- pd |> 
    select(cell_type_class, cell_type_broad, cell_type_fine, passALL_metricQC, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(cell_type_fine,cell_type_broad, cell_type_class,passALL_metricQC), names_to = "metric") |>
    mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))) |>
    ggplot() +
    geom_violin(aes(x = cell_type_fine, y = value, fill = cell_type_fine, colour = passALL_metricQC), 
                scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    geom_hline(data = cell_class_cutoffs, aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    geom_label(data = cell_class_cutoffs, aes(y = cutoff, label = cutoff_anno), x = 3, color = "blue", vjust = -.5) +
    facet_grid(metric~cell_type_class, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Quality Metric Value", x = "preliminary cell type cluster, snn k=10")

ggsave(qc_violin_plot_all, filename = here(plot_dir, "sn_clusterQC-QCmetricViolin_ALL.png"), width = 12, height = 12)


#### n cell barplot ####
n_cell_barplot <- cluster_info |>
    ggplot(aes(x = cell_type_fine, y = n, fill = cell_type_fine, color = passALL_metricQC)) +
    geom_col() +
    geom_text(aes(label = n), size = 2, vjust = -0.5) +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(n_cell_barplot, filename = here(plot_dir, "sn_clusterQC-barplot_ct_n.png"), width = 10)

####  Cell Type Proportions ####

cell_class_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, cell_type_class) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n))

sample_neuron_order <- cell_class_proportions |> filter(cell_type_class == "neuron") |> arrange(prop) |> pull(sample_id)

cell_type_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, cell_type_fine) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = sample_neuron_order),
           passALL_metricQC = TRUE) #TODO fix passALL_metricQC bug
    # ungroup() |>
    # left_join(cluster_metrics_pass, by = join_by(cell_type_fine))

cell_type_proportions |> filter(!passALL_metricQC) |> arrange(-prop)

cell_type_proportion_bar <- cell_type_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = cell_type_fine)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), ""),
                  color = passALL_metricQC),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red"), guide = "none") +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(cell_type_proportion_bar, filename = here(plot_dir, "sn_clusterQC-barplot_ct_prop.png"), width = 10)

## evaluate proportion sample by cell 

## create sample color pallet
load(here(here("processed-data", "project_colors.Rdata")), verbose = TRUE)
apoe_colors_tb <- pd |> 
    select(sample_id, APOE) |> 
    unique() |>
    group_by(APOE) |>
    arrange(sample_id) |>
    mutate(APOE_color = paste0(APOE, ".", row_number())) |>
    arrange(APOE_color)

apoe_colors <- DeconvoBuddies::create_cell_colors(apoe_colors_tb$APOE_color, pallet = APOE_genotype_colors, split = "\\.")
sample_colors <- apoe_colors$fine
names(sample_colors) <- apoe_colors_tb$sample_id

sample_proportions <- pd |>
    group_by(sample_id, cell_type_fine) |>
    summarize(n = n()) |> 
    group_by(cell_type_fine) |>
    mutate(prop = n/sum(n),
          sample_id = factor(sample_id, levels = apoe_colors_tb$sample_id)) |>
    left_join(cluster_metrics_pass)

## check for cell types from few samples 
sample_proportions |> filter(prop > 0.3)

sample_qc <- sample_proportions |> 
    group_by(cell_type_fine, passALL_metricQC) |>
    summarize(n_samples = length(unique(sample_id)),
              max_prop = max(prop)) |>
    mutate(pass_sampleQC = max_prop < 0.5 & n_samples > 5)

sample_qc |> arrange(n_samples)
sample_qc |> arrange(-max_prop)

## plot scatter
n_samples_vs_max_prop <- sample_qc |>
    ggplot(aes(n_samples, max_prop, fill = cell_type_fine, shape = passALL_metricQC, color = pass_sampleQC)) +
    geom_point() +
    geom_text_repel(aes(label = cell_type_fine), size = 2) +
    geom_hline(yintercept = 0.5, color = "blue", linetype = "dashed") +
    geom_vline(xintercept = 5, color = "blue", linetype = "dashed") +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24)) +
    theme_bw()

ggsave(n_samples_vs_max_prop, filename = here(plot_dir, "sn_clusterQC-scater_n_samples_vs_max_prop.png"))

## plot sample by cell type proportions bar
sample_proportion_bar <- sample_proportions |>
    ggplot(aes(x = cell_type_fine, y = prop, fill = sample_id)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > 0.3, sprintf("%s - %g", sample_id, round(prop, 3)), "")),
              position = position_stack(vjust = .5),
              size = 2, angle = 90) +
    theme_bw() +
    scale_fill_manual(values = sample_colors) +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(~passALL_metricQC,  scales = "free_x", space = "free_x")

ggsave(sample_proportion_bar, filename = here(plot_dir, "sn_clusterQC-barplot_sample_prop.png"), width = 10)

sample_n_bar <- sample_proportions |>
    ggplot(aes(x = cell_type_fine, y = n, fill = sample_id)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > 0.3, sprintf("%s - %g", sample_id, n), "")),
              position = position_stack(vjust = .5),
              size = 2.3, angle = 90) +
    theme_bw() +
    scale_fill_manual(values = sample_colors) +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(~passALL_metricQC,  scales = "free_x", space = "free_x")

ggsave(sample_n_bar, filename = here(plot_dir, "sn_clusterQC-barplot_sample_n.png"), width = 10)

#### Add sample fraction QC to cluster info and save ####

cluster_info <- cluster_info |>
    left_join(sample_qc) |>
    mutate(pass_clusterQC = passALL_metricQC & pass_sampleQC)

cluster_info |> ungroup()|> count(pass_clusterQC, passALL_metricQC, pass_sampleQC)

cluster_info |>
    group_by(pass_clusterQC) |>
    summarize(n = sum(n),
              prop = n/nrow(pd),
              n_clusters = length(unique(cell_type_fine)))


write.csv(cluster_info, file = here(data_dir, "sn_clusterQC_info.csv"), row.names = FALSE)

#### plot reduced dims ####

pd <- pd |> left_join(sample_qc) |> mutate(pass_clusterQC = passALL_metricQC & pass_sampleQC)
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
                                              my_var = "pass_clusterQC", 
                                              color_pal = c(`TRUE` = "black", `FALSE` = "red"),
                                              sufix = "prelim_k10"))

## filtered 
sce <- sce[,sce$pass_clusterQC]

walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn_clusterQCpass",
                                              var_type = "cat", 
                                              dimred = .x, 
                                              my_var = "cell_type_broad", 
                                              color_pal = cell_type_colors$broad,
                                              sufix = "prelim_k10"))


# slurmjobs::job_single('09_cluster_QC',
#                       create_shell = TRUE, memory = '10G',
#                       command = "Rscript 09_cluster_QC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

