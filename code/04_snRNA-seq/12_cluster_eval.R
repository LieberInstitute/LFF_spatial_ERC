## Louise Huuki-Myers, February 2024
## Evaluate different clustering results with Silhouette width and Jacquard matrix

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

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "09_cluster_eval")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "09_cluster_eval")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_harmony"))
sce

#### load cluster outputs ####
message(Sys.time(), " - Load cluster data")
cluster_fn <- list.files(here("processed-data", "04_snRNA-seq", "07_cluster_sn"), full.names = TRUE)  
names(cluster_fn) <- ss(basename(cluster_fn), "_", 3)

clusters <- map(cluster_fn, ~get(load(.x)))

map_int(clusters, max)
# k10 k15 k20 
# 62  28  33 

map_int(clusters, ~sum(table(.x) < 100))
# k10 k15 k20 
# 12   4   1 

map(clusters, table)

sce$snn_k10 <- sprintf("k10c%02d", clusters$k10)
sce$snn_k15 <- sprintf("k15c%02d", clusters$k15)
sce$snn_k20 <- sprintf("k20c%02d", clusters$k20)

#### load sc-type output ####
message(Sys.time(), " - scType data")

sctype_fn <- list.files(here("processed-data", "04_snRNA-seq", "10_sctype"), pattern = ".csv", full.names = TRUE)
names(sctype_fn) <- ss(basename(sctype_fn), "_|\\.", 4)

sctype <- map(sctype_fn, ~read.csv(.x, row.names = 1))

map(sctype, ~.x |> group_by(type) |>  summarize(ncells = sum(ncells), n_cluster = n()))
# $k15
# # A tibble: 7 × 3
# type                            ncells n_cluster
# <chr>                            <int>     <int>
# 1 Astrocytes                       24387         6
# 2 Endothelial cells                 1739         1
# 3 GABAergic neurons                11134         3
# 4 Glutamatergic neurons            19241        11
# 5 Microglial cells                 15441         2
# 6 Oligodendrocyte precursor cells  13070         1
# 7 Oligodendrocytes                 55107         4

sort(unique(unlist(map(sctype, "type"))))

## short names 
type_short <- tibble(type = c(
    'Astrocytes',
    'Endothelial cells',
    'Microglial cells',
    'Oligodendrocytes',
    'Oligodendrocyte precursor cells',
    'Glutamatergic neurons',
    'GABAergic neurons',
    'Mature neurons'), 
    type_short  = c(
        'Astro',
        'Endo',
        "Micro",
        'Oligo',
        'OPC',
        'Excit',
        'Inhib',
        'Neu')
)

cell_type_levels = type_short$type_short

sctype <- map(sctype , ~.x|> 
                  left_join(type_short) |>
                  group_by(type) |>
                  arrange(-ncells) |> 
                  mutate(type_rank = row_number(),
                         ct_fine = paste0(type_short, 
                                          ".",
                                          width = str_pad(type_rank, 
                                                  nchar(as.character(max(type_rank))),
                                          side = "left",
                                          pad = "0")),
                         type_short = factor(type_short, levels = cell_type_levels),
                         )
              )

map(sctype, ~levels(.x$type_short))

fine_levels <- map(sctype, ~.x |> arrange(type_short) |> pull(ct_fine))

sctype <- map2(sctype, fine_levels, ~.x |>
                   mutate(ct_fine = factor(ct_fine, levels = .y)) |>
                   arrange(ct_fine)
)

map(sctype, ~levels(.x$ct_fine))

levels(sctype$k15$type_short)

map(sctype, ~.x |> arrange(-type_rank))
map(sctype, ~.x |> count(type_short) |> arrange(type_short))

# $k15
# # A tibble: 7 × 3
# # Groups:   type [7]
# type                            type_short     n
# <chr>                           <fct>      <int>
# 1 Astrocytes                      Astro          6
# 2 Endothelial cells               Endo           1
# 3 Microglial cells                Micro          2
# 4 Oligodendrocytes                Oligo          4
# 5 Oligodendrocyte precursor cells OPC            1
# 6 GABAergic neurons               Neu_GABA       3
# 7 Glutamatergic neurons           Neu_Glut      11

## write annotation csv 
walk2(sctype, names(sctype), ~write.csv(.x, file = here(data_dir, paste0("cluster_annotation_snn_", .y,".csv")), row.names = FALSE))

## match to clusters in sce
anno_table <- map2_dfc(sctype, names(sctype), function(sct, name){
    cl_col <- paste0("snn_", name)
    anno_df = data.frame(
        ct_broad = factor(sct$type_short[match(sce[[cl_col]], sct$cluster)],
                          cell_type_levels),
        ct_fine = factor(sct$ct_fine[match(sce[[cl_col]], sct$cluster)], 
                         fine_levels[[name]])
    )
    
    colnames(anno_df) = paste0(colnames(anno_df), "_", name)
    
    return(anno_df)
})

head(anno_table)

table(anno_table$ct_broad_k20, anno_table$ct_broad_k15)

levels(anno_table$ct_fine_k15)
levels(anno_table$ct_fine_k20)

## add to colData
colData(sce) <- cbind(colData(sce), anno_table)

#### Define colors for fine cell types ####

ks <- c("k10", "k15", "k20")
names(ks) <- ks

cell_type_colors <- map(ks ,~create_cell_colors(cell_types = levels(anno_table[[paste0("ct_fine_", .x)]]),
                                                pallet_name = "classic", 
                                                split = "\\."))
# $k20
# $broad
# Astro      Endo     Micro     Oligo       OPC     Excit     Inhib 
# "#3BB273" "#FF56AF" "#663894" "#F57A00" "#D2B037" "#247FBC" "#E83E38" 
# 
# $fine
# Astro.1   Astro.2   Astro.3   Astro.4   Astro.5    Endo.1   Micro.1   Micro.2   Micro.3   Micro.4   Oligo.1   Oligo.2 
# "#3BB273" "#62C18F" "#89D0AB" "#B0E0C7" "#D7EFE3" "#FF56AF" "#663894" "#8C69AE" "#B29BC9" "#D8CDE4" "#F57A00" "#F79433" 
# Oligo.3   Oligo.4   Oligo.5     OPC.1  Excit.01  Excit.02  Excit.03  Excit.04  Excit.05  Excit.06  Excit.07  Excit.08 
# "#F9AF66" "#FBC999" "#FDE4CC" "#D2B037" "#247FBC" "#3388C0" "#4391C5" "#529ACA" "#62A3CF" "#72ACD3" "#81B5D8" "#91BFDD" 
# Excit.09  Excit.10  Excit.11  Excit.12  Excit.13  Excit.14   Inhib.1   Inhib.2   Inhib.3 
# "#A1C8E2" "#B0D1E7" "#C0DAEB" "#D0E3F0" "#DFECF5" "#EFF5FA" "#E83E38" "#EF7E7A" "#F7BEBC"

## save all color output
save(cell_type_colors, file = here(data_dir, "cell_type_colors_allK.Rdata"))

#### plot on UMAP + TSNE ####
message(Sys.time(), " - Plot UMAP + TSNE")

## plot clusters
walk(names(clusters), function(k){
    walk(c("UMAP", "TSNE"),  
         ~my_plot_reduced_dim(sce, 
                              prefix = "sn", 
                              var_type = "cat", 
                              dimred = .x, 
                              my_var = paste0("snn_",k), 
                              sufix = "HARMONY"))
})
## plot annotations
walk2(colnames(anno_table), unlist(cell_type_colors, recursive = FALSE), function(ct, colors){
    
    walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn",
                                                  var_type = "cat", 
                                                  dimred = .x, 
                                                  my_var = ct, 
                                                  color_pal = colors,
                                                  sufix = "HARMONY"))
})


table(sce$snn_k15, sce$snn_k20)

#### plot marker genes ####
message(Sys.time(), " - Plot marker genes")
## summarize marker genes from lit 
lit_markers <- read.csv(here("processed-data", "04_snRNA-seq", "lit_marker_genes.csv")) |>
    group_by(gene_name, cell_type) |>
    summarize(n_studies = n(),
            studies = paste0(source, collapse = ",")) |>
    mutate(in_data = gene_name %in% rowData(sce)$Symbol)
    
lit_markers |> write_csv(here("processed-data", "04_snRNA-seq", "lit_marker_summary.csv"))

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
                         pdf_fn = here(plot_dir, "ssn_k15_lit_markers.pdf"),
                         cellType_col = "ct_fine_k15",
                         gene_name_col = "Symbol",
                         color_pal = cell_type_colors$k15$fine,
                         )

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ssn_k20_lit_markers.pdf"),
                         cellType_col = "ct_fine_k20",
                         gene_name_col = "Symbol",
                         color_pal = cell_type_colors$k20$fine,
                         )

#### Compare clustering w/ Jaccard Index ####
message(Sys.time(), " - Compare clusters")

walk(c("broad", "fine"), function(resolution){
    
    cluster_names <- sprintf("ct_%s_k%d", resolution, c(10, 15, 20))
    
    cluster_combos <- as.data.frame(t(combn(cluster_names, 2)))
    colnames(cluster_combos) <- c("c1", "c2")
    
    cluster_combos <- cluster_combos |> 
        mutate(name = gsub(sprintf("ct_%s_", resolution), "",paste0(c1, "_",c2))) 
    
    jacc.mat <- map2(cluster_combos$c1, cluster_combos$c2, ~linkClustersMatrix(sce[[.x]], sce[[.y]]))
    names(jacc.mat) <- cluster_combos$name
    
    walk2(jacc.mat, names(jacc.mat), function(j, name){
        pdf(here(plot_dir, sprintf("jacc_matrix_ct_%s_%s.pdf", resolution, name)), height = 12, width = 8)
        print(Heatmap(j,
                      name = "Correspondence",
                      col = c("black", viridisLite::plasma(100)),
                      na_col = "black"
        ))
        dev.off()
    })
    
})

## examine hierarchical clustering ####
message(Sys.time(), " - Plot HC with ct names")

h_clus_fn <- list.files(here("processed-data", "04_snRNA-seq", "08_hierarchical_cluster"), full.names = TRUE)
names(h_clus_fn) <- ss(basename(h_clus_fn), "_|\\.", 3)

h_clus <- map(h_clus_fn, ~mget(load(.x)))

dend_ct <- map2(h_clus, sctype, function(hc, sct){
    
    d <- hc$dend
    labels(d) <- sct$ct_fine[match(labels(d), sct$cluster)]
    return(d)
})

labels(dend_ct$k10)

walk2(dend_ct, names(dend_ct), function(d, name){
    pdf(here(plot_dir, sprintf("dend_ct_%s.pdf", name)), height = 10, width = 12)
    par(cex = 0.6, font = 2)
    plot(d, main = sprintf("hierarchical cluster dend - %s", name), horiz = TRUE)
    dev.off()
    })


#### Add colors to metadata ####

metadata(sce)$cell_type_colors <- cell_type_colors$k20

#### Output annotations ####
message(Sys.time(), " - Save annotated sce")

save(sce, file = here("processed-data", "spe_objects", "sce_ERC.Rdata"))
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce_objects", "sce_ERC"), replace=TRUE)

# slurmjobs::job_single('09_cluster_eval', create_shell = TRUE, memory = '5G', command = "09_cluster_eval.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

