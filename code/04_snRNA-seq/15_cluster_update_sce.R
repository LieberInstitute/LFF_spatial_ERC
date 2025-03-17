## Louise Huuki-Myers, March 2025
## Add cluster data to sce object, explore and annotate clusters

library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
library("ComplexHeatmap")
library("bluster")
library("dendextend")
library("readxl")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "15_cluster_update_sce")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "15_cluster_update_sce")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_reprocess"))
sce

## clean up colData
colData(sce) <- colData(sce)[,!grepl("snn|ct_", colnames(colData(sce)))]

#### Add optimal cluster output ####
optimal_k = 13
message("Load optimal clusters from K=", optimal_k)
load(here("processed-data", "04_snRNA-seq", "11_cluster_sn", sprintf("walktrap_snn_k%d_clusters.Rdata", optimal_k)), verbose = TRUE)
# clusters

sce$snn_kOpt <- sprintf("k%dc%02d", optimal_k, clusters)

#### load sc-type output ####
load(here("processed-data", "04_snRNA-seq", "13_sctype_final", "sctype_final-sctype.Rdata"), verbose = TRUE)

anno_notes <- readxl::read_excel(here("processed-data", "04_snRNA-seq", "13_sctype_final", "sctype_final-NOTES.xlsx"))

#sctype
table(sctype$cell_type_broad)

anno_table <- sctype |> 
    ungroup() |>
    left_join(anno_notes |> select(cluster, cell_type_anno = guess)) |>
    select(cluster, cell_type_class, cell_type_broad, cell_type_fine, cell_type_anno) |>
    column_to_rownames("cluster")

cell_type_anno_levels <- anno_table |> arrange(cell_type_broad, cell_type_anno) |> pull(cell_type_anno)

anno_table$cell_type_anno <- factor(anno_table$cell_type_anno, levels = cell_type_anno_levels)

levels(anno_table$cell_type_broad)
levels(anno_table$cell_type_fine)
levels(anno_table$cell_type_anno)

anno_table <- anno_table[sce$snn_kOpt,]
dim(anno_table)

#### annotation notes ####

## add to colData
colData(sce) <- cbind(colData(sce), anno_table)

#### Define colors for fine cell types ####
cell_type_colors <- DeconvoBuddies::create_cell_colors(cell_types = levels(sctype$cell_type_fine),
                                                pallet_name = "classic", 
                                                split = "\\.")
# $broad
# Astro      Endo     Micro     Oligo       OPC     Excit     Inhib       Neu 
# "#3BB273" "#FF56AF" "#663894" "#F57A00" "#D2B037" "#247FBC" "#E83E38" "#4E586A" 
# 
# $fine
# Astro.1   Astro.2    Endo.1   Micro.1   Micro.2   Oligo.1   Oligo.2   Oligo.3     OPC.1   Excit.1   Excit.2   Excit.3 
# "#3BB273" "#9DD8B9" "#FF56AF" "#663894" "#B29BC9" "#F57A00" "#F8A655" "#FBD2AA" "#D2B037" "#247FBC" "#3C8DC3" "#549BCA" 
# Excit.4   Excit.5   Excit.6   Excit.7   Excit.8   Excit.9   Inhib.1   Inhib.2   Inhib.3     Neu.1     Neu.2 
# "#6DA9D2" "#85B7D9" "#9DC6E1" "#B5D4E8" "#CEE2F0" "#E6F0F7" "#E83E38" "#EF7E7A" "#F7BEBC" "#4E586A" "#A6ABB4" 

cell_type_colors_anno <- DeconvoBuddies::create_cell_colors(cell_types = cell_type_anno_levels,
                                                       pallet_name = "classic", 
                                                       split = "\\.")
# $broad
# Astro      Endo     Micro       OPC     Oligo     Excit     Inhib 
# "#3BB273" "#FF56AF" "#663894" "#F57A00" "#D2B037" "#247FBC" "#E83E38" 
# 
# $fine
# Astro.1          Astro.2             Endo          Micro.1          Micro.2          Oligo.1          Oligo.2 
# "#3BB273"        "#9DD8B9"        "#FF56AF"        "#663894"        "#B29BC9"        "#F57A00"        "#F8A655" 
# Oligo.3              OPC     Excit.L2_3.1     Excit.L2_3.2     Excit.L2_3.3       Excit.L4_5         Excit.L5 
# "#FBD2AA"        "#D2B037"        "#247FBC"        "#3F8FC4"        "#5A9FCC"        "#76AFD5"        "#91BFDD" 
# Excit.L5_6_NP      Excit.L6_CT        Excit.L6b       Inhib.Pax6 Inhib.Lamp5_Lhx6      Inhib.Pvalb        Inhib.Vip 
# "#ACCFE5"        "#C8DFEE"        "#E3EFF6"        "#E83E38"        "#EB5E59"        "#EF7E7A"        "#F39E9B" 
# Inhib.Chandelier        Inhib.Sst 
# "#F7BEBC"        "#FBDEDD" 

cell_type_colors$anno <- cell_type_colors_anno$fine

## save all color output
save(cell_type_colors, file = here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"))

#### plot on UMAP + TSNE ####
message(Sys.time(), " - Plot UMAP + TSNE")

## plot annotated clusters
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(sce,
                          prefix = "ERC_sn_reprocess",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "cell_type_fine",
                          color_pal = cell_type_colors$fine,
                          suffix = "sctype"))

walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(sce,
                          prefix = "ERC_sn_reprocess",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "cell_type_fine",
                          color_pal = cell_type_colors$anno,
                          suffix = "sctype"))


#### plot marker genes ####
message(Sys.time(), " - Plot marker genes")

## read in marker genes from lit
lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    mutate(in_data = gene_name %in% rowData(sce)$Symbol,
           cell_type_broad = factor(cell_type_broad, levels = levels(anno_table$cell_type_broad))) |>
    arrange(cell_type_broad)

## missing from our data
lit_markers |> filter(!in_data)
# gene_name cell_type          n_studies studies                  cell_type_broad in_data
# <chr>     <chr>                  <dbl> <chr>                    <chr>           <lgl>  
# 1 CEMP      Fibroblast                 1 Davila-Velderrain et al. Endo            FALSE  
# 2 GR1A1     Glutamate receptor         1 Grubman et al.           Excit           FALSE  
# 3 PDCH15    OPC                        1 Grubman et al.           OPC             FALSE 

lit_markers <- lit_markers |> filter(in_data)

lit_markers_list <- map(splitit(lit_markers$cell_type), ~lit_markers$gene_name[.x])

lit_marker_order <- lit_markers |> count(cell_type_broad, cell_type) |> pull(cell_type)
lit_markers_list <- lit_markers_list[lit_marker_order]

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_sn_sctype_lit_markers.pdf"),
                         cellType_col = "cell_type_anno",
                         gene_name_col = "Symbol",
                         color_pal = cell_type_colors$anno,
                         )

#### summary by cluster ####
message(Sys.time(), " - Summarize cluster")

pd <- as.data.frame(colData(sce))

cluster_info <- pd |>
    group_by(snn_kOpt, cell_type_class, cell_type_broad, cell_type_fine, cell_type_anno) |>
    summarize(n = n(),
              prop = n/ncol(sce),
              median_sum = median(sum),
              median_detected = median(detected),
              median_Mito_percent = median(subsets_Mito_percent),
              median_scDblFinder.score = median(scDblFinder.score)
    )

cluster_info |> arrange(median_sum)

cell_class_cutoffs <- read.csv(here("processed-data", "04_snRNA-seq", "09_cluster_QC","ERC_sn_cell_class_cutoffs.csv"))

cluster_metrics_long <- pd |> 
    select(snn_kOpt, cell_type_class, cell_type_broad, cell_type_fine, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(snn_kOpt, cell_type_fine,cell_type_broad, cell_type_class), names_to = "metric") |>
    group_by(snn_kOpt, cell_type_fine, cell_type_class, cell_type_broad, metric) |>
    summarize(median = median(value)) |>
    left_join(cell_class_cutoffs |> select(cell_type_class, metric, cutoff, cutoff_anno)) |>
    mutate(pass_metricQC = case_when(metric %in% c("scDblFinder.score", "subsets_Mito_percent") ~ median < cutoff,
                                     TRUE ~ median > cutoff),
           metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))
    )|>
    group_by(snn_kOpt, cell_type_fine, cell_type_class, cell_type_broad) |>
    mutate(passALL_metricQC = all(pass_metricQC))

cluster_metrics_long |> ungroup() |> count(pass_metricQC, passALL_metricQC)

cluster_metrics_pass <- cluster_metrics_long |> ungroup() |> select(cell_type_fine, passALL_metricQC)  |> unique()
cluster_metrics_pass |> filter(!passALL_metricQC)

## add to info
pd <- pd |> left_join(cluster_metrics_pass)
cluster_info <- cluster_info |> left_join(cluster_metrics_pass)

write.csv(cluster_info, file = here(data_dir, "ERC_sn_cluster_info.csv"), row.names = FALSE)

#### Plot Quality Metrics ####

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

ggsave(qc_median_boxplot, filename = here(plot_dir, "ERC_sn_QCmedian_boxplot.png"))    

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

ggsave(qc_violin_plot_all, filename = here(plot_dir, "ERC_sn_QCmetricViolin_ALL.png"), width = 12, height = 12)

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

ggsave(n_cell_barplot, filename = here(plot_dir, "ERC_sn_barplot_ct_n_qc.png"), width = 10)

## with out QC fail clusters
n_cell_barplot <- cluster_info |>
    filter(passALL_metricQC) |>
    ggplot(aes(x = cell_type_fine, y = n, fill = cell_type_fine, color = passALL_metricQC)) +
    geom_col() +
    geom_text(aes(label = n), size = 2, vjust = -0.5) +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(n_cell_barplot, filename = here(plot_dir, "ERC_sn_barplot_ct_n.png"), width = 10)

#### Drop problamatic cluster ####
colData(sce) <- DataFrame(pd)
table(sce$passALL_metricQC)

sce <- sce[,sce$passALL_metricQC]
dim(sce)

sce$cell_type_fine <- droplevels(sce$cell_type_fine)

levels(sce$cell_type_fine)
pd <- as.data.frame(colData(sce))
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
           sample_id = factor(sample_id, levels = sample_neuron_order))

cell_type_proportion_bar <- cell_type_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = cell_type_fine)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors$fine, guide = "none") +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(cell_type_proportion_bar, filename = here(plot_dir, "ERC_sn_barplot_ct_prop.png"), width = 10)

## evaluate proportion sample by cell 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
# we need sample_colors

sample_proportions <- pd |>
    group_by(sample_id, cell_type_fine) |>
    summarize(n = n()) |> 
    group_by(cell_type_fine) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = names(sample_colors))) 

## check for cell types from few samples 
sample_proportions |> arrange(-prop)

sample_qc <- sample_proportions |> 
    group_by(cell_type_fine) |>
    summarize(n_samples = length(unique(sample_id)),
              max_prop = max(prop)) |>
    mutate(pass_sampleQC = max_prop < 0.5 & n_samples > 5)

## lowest contibution is 19 samples 
sample_qc |> arrange(n_samples)
sample_qc |> arrange(-max_prop)

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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sample_proportion_bar, filename = here(plot_dir, "ERC_sn_barplot_sample_prop.png"), width = 10)

    
# #### Compare clustering w/ Jaccard Index ####
# message(Sys.time(), " - Compare clusters")
# 
# walk(c("broad", "fine"), function(resolution){
#     
#     cluster_names <- sprintf("ct_%s_k%d", resolution, c(10, 15, 20))
#     
#     cluster_combos <- as.data.frame(t(combn(cluster_names, 2)))
#     colnames(cluster_combos) <- c("c1", "c2")
#     
#     cluster_combos <- cluster_combos |> 
#         mutate(name = gsub(sprintf("ct_%s_", resolution), "",paste0(c1, "_",c2))) 
#     
#     jacc.mat <- map2(cluster_combos$c1, cluster_combos$c2, ~linkClustersMatrix(sce[[.x]], sce[[.y]]))
#     names(jacc.mat) <- cluster_combos$name
#     
#     walk2(jacc.mat, names(jacc.mat), function(j, name){
#         pdf(here(plot_dir, sprintf("jacc_matrix_ct_%s_%s.pdf", resolution, name)), height = 12, width = 8)
#         print(Heatmap(j,
#                       name = "Correspondence",
#                       col = c("black", viridisLite::plasma(100)),
#                       na_col = "black"
#         ))
#         dev.off()
#     })
#     
# })

# ## examine hierarchical clustering ####
# message(Sys.time(), " - Plot HC with ct names")
# 
# h_clus_fn <- list.files(here("processed-data", "04_snRNA-seq", "08_hierarchical_cluster"), full.names = TRUE)
# names(h_clus_fn) <- ss(basename(h_clus_fn), "_|\\.", 3)
# 
# h_clus <- map(h_clus_fn, ~mget(load(.x)))
# 
# dend_ct <- map2(h_clus, sctype, function(hc, sct){
#     
#     d <- hc$dend
#     labels(d) <- sct$ct_fine[match(labels(d), sct$cluster)]
#     return(d)
# })
# 
# labels(dend_ct$k10)
# 
# walk2(dend_ct, names(dend_ct), function(d, name){
#     pdf(here(plot_dir, sprintf("dend_ct_%s.pdf", name)), height = 10, width = 12)
#     par(cex = 0.6, font = 2)
#     plot(d, main = sprintf("hierarchical cluster dend - %s", name), horiz = TRUE)
#     dev.off()
#     })


#### Add colors to metadata ####

metadata(sce)$cell_type_colors <- cell_type_colors

#### Output annotations ####
message(Sys.time(), " - Save annotated sce")

# save(sce, file = here("processed-data", "spe_objects", "sce_ERC.Rdata"))
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce_objects", "sce_ERC"), replace=TRUE)

# slurmjobs::job_single('15_cluster_update_sce', create_shell = TRUE, memory = '10G', command = "15_cluster_update_sce.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

