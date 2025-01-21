## Louise Huuki-Myers, September 2024
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

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "09_cluster_annotation")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "09_cluster_annotation")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_harmony"))
sce

## load cluster outputs
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
    'GABAergic neurons',
    'Glutamatergic neurons',
    'Mature neurons'), 
    type_short  = c(
        'Astro',
        'Endo',
        "Micro",
        'Oligo',
        'OPC',
        'Neu_GABA',
        'Neu_Glut',
        'Neu_matu')
)

cell_type_levels = type_short$type_short

sctype <- map(sctype , ~.x|> 
                  left_join(type_short) |>
                  group_by(type) |>
                  arrange(-ncells) |> 
                  mutate(type_rank = row_number(),
                         ct_fine = paste0(type_short, 
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
         mutate(ct_fine = factor(ct_fine, levels = .y))
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
        ct_fine = factor(sct$ct_fine[match(sce[[cl_col]], sct$cluster)]), 
        fine_levels[[.y]])
    
    colnames(anno_df) = paste0(colnames(anno_df), "_", name)
    
    return(anno_df)
})

head(anno_table)

table(anno_table$ct_broad_k20, anno_table$ct_broad_k15)

levels(anno_table$ct_fine_k15)

## add to colData
colData(sce) <- cbind(colData(sce), anno_table)

## plot on UMAP + TSNE
walk(names(clusters), function(k){
    walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn", var_type = "cat", dimred = .x, my_var = paste0("snn_",k), sufix = "HARMONY"))
})

walk(colnames(anno_table), function(ct){
    walk(c("UMAP", "TSNE"),  ~my_plot_reduced_dim(sce, prefix = "sn", var_type = "cat", dimred = .x, my_var = ct, sufix = "HARMONY"))
})


table(sce$snn_k15, sce$snn_k20)

#### plot marker genes ####

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
                         gene_name_col = "Symbol"
                         )

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ssn_k20_lit_markers.pdf"),
                         cellType_col = "ct_fine_k20",
                         gene_name_col = "Symbol"
                         )



#### Compare clustering w/ Jaccard Index ####

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

## examine heriachial clustering ####

# slurmjobs::job_single('09_cluster_annotation', create_shell = TRUE, memory = '5G', command = "09_cluster_annotation.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

