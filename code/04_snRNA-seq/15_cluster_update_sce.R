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
#sctype

anno_table <- sctype |> 
    ungroup() |>
    select(cluster, cell_type_class, cell_type_broad, cell_type_fine) |>
    column_to_rownames("cluster")

levels(anno_table$cell_type_broad)
levels(anno_table$cell_type_fine)

anno_table <- anno_table[sce$snn_kOpt,]
dim(anno_table)

## add to colData
colData(sce) <- cbind(colData(sce), anno_table)

#### Define colors for fine cell types ####
cell_type_colors <- create_cell_colors(cell_types = levels(sctype$cell_type_fine),
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


## save all color output
save(cell_type_colors, file = here(data_dir, "cell_type_colors.Rdata"))

#### plot on UMAP + TSNE ####
message(Sys.time(), " - Plot UMAP + TSNE")

## plot clusters
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(sce,
                          prefix = "ERC_sn_reprocess",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "cell_type_fine",
                          color_pal = "cell_type_colors"
                          suffix = "sctype"))

## plot annotations
walk2(c("cell_type_broad", "cell_type_fine"), cell_type_colors, function(ct, colors){
    
    walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "ERC_sn_reprocess",
                                                  var_type = "cat", 
                                                  dimred = .x, 
                                                  my_var = ct, 
                                                  color_pal = colors,
                                                  suffix = "sctype"))
})


#### plot marker genes ####
message(Sys.time(), " - Plot marker genes")

## read in marker genes from lit
lit_markers <-  read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv"))

## missing from our data
lit_markers |> filter(!in_data)
# gene_name cell_type          n_studies studies        in_data
# <chr>     <chr>                  <int> <chr>          <lgl>  
# 1 GR1A1     Glutamate receptor         1 Grubman et al. FALSE  
# 2 PDCH15    OPC                        1 Grubman et al. FALSE

lit_markers <- lit_markers |> filter(in_data)

lit_markers_list <- map(splitit(lit_markers$cell_type), ~lit_markers$gene_name[.x])

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_sn_sctype_lit_markers.pdf"),
                         cellType_col = "cell_type_fine",
                         gene_name_col = "Symbol",
                         color_pal = cell_type_colors$fine,
                         )

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

save(sce, file = here("processed-data", "spe_objects", "sce_ERC.Rdata"))
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce_objects", "sce_ERC"), replace=TRUE)

# slurmjobs::job_single('15_cluster_update_sce', create_shell = TRUE, memory = '5G', command = "15_cluster_update_sce.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

