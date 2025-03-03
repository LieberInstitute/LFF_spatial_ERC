## Louise Huuki-Myers, January 2025
## Annotate preliminary clusters with SC Type

## load libraries
library("tidyverse")
library("Seurat")
library("HDF5Array")
library("HGNChelper")
library("openxlsx")
library("here")
library("sessioninfo")
library("bluster")

## source sc-type functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")

## Prep directories
# plot_dir <- here("plots", "04_snRNA-seq", "13_sctype_final")
# if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "13_sctype_final")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load sce object
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_harmony"))
sce

## load cluster outputs
load(here("processed-data", "04_snRNA-seq", "07_cluster_sn_prelim", "walktrap_snn_k10_clusters_prelim.Rdata"), verbose = TRUE)
# clusters

sce$snn_k10 <- sprintf("k10c%02d", clusters)

message("k10 prelim clusters: ", length(unique(sce$snn_k10)))

table(sce$snn_k10)


# table(sce$quick_cluster, sce$snn_k10)

## Adjusted Rand Index
# 0.5 corresponds to “good” similarity
message("Rand Index w/ Quick Cluster")
pairwiseRand(sce$quick_cluster, sce$snn_k10, mode = "index") #0.4513222


#### SC Type #### 
# get cell-type-specific gene sets from our in-built database (DB)
# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# prepare gene sets
gs_list <- gene_sets_prepare(db_, "Brain")

rownames(sce) <- rowData(sce)$Symbol

message(Sys.time(), "- sctype score")
es.max <- sctype_score(scRNAseqData = as.matrix(logcounts(sce)), 
                       scaled = TRUE, 
                       gs = gs_list$gs_positive)

dim(es.max)

es.max[1:5, 1:5]
#                          cell1     cell2     cell3    cell4     cell5
# Astrocytes           1.0788117 0.7292370 0.7731955 6.940781 0.3259594
# Cholinergic neurons  0.0000000 0.2280047 0.0000000 0.000000 0.0000000
# Dopaminergic neurons 0.0000000 0.8389230 1.0932569 0.000000 0.2723543
# Endothelial cells    0.7930059 0.0000000 0.0000000 0.000000 0.0000000
# GABAergic neurons    1.4841430 2.3017959 4.7711181 1.622680 5.6416159

save(es.max, file = here(data_dir, "sctype_es_max.Rdata"))

#### annotate clusters ####
message(Sys.time(), " - Annotate clusters")

cl_col <- "snn_k10"

clusters <- sort(unique(sce[[cl_col]]))

## compile scores
cL_results <- purrr::map_dfr(clusters, function(cluster){
    cluster_index = sce[[cl_col]] == cluster
    
    es.max.cl = sort(rowSums(es.max[,cluster_index]), decreasing = !0)
    cL_resutls <- head(tibble(cluster = cluster, 
                              type = names(es.max.cl), 
                              scores = es.max.cl, 
                              ncells = sum(cluster_index)
    ),10)
    return(cL_resutls)
})

sctype_scores <-  cL_results |> 
    group_by(cluster) |> 
    top_n(n = 1, wt = scores)  |>
    mutate(confident = scores >= (ncells/4))

#### Annotate ####

sctype_scores |> group_by(type) |>  summarize(ncells = sum(ncells), n_cluster = n())
# type                            ncells n_cluster
# <chr>                            <int>     <int>
# 1 Astrocytes                       24913        11
# 2 Endothelial cells                 1693         1
# 3 GABAergic neurons                10668         6
# 4 Glutamatergic neurons            19327        23
# 5 Mature neurons                     404         1
# 6 Microglial cells                 15223         6
# 7 Oligodendrocyte precursor cells  15373         5
# 8 Oligodendrocytes                 52518         9

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
    cell_type_broad  = c(
        'Astro',
        'Endo',
        "Micro",
        'Oligo',
        'OPC',
        'Excit',
        'Inhib',
        'Neu'),
    cell_type_class =c(
        rep("glia", 5),
        rep("neuron", 3)
    )
)


cell_type_levels = type_short$cell_type_broad

sctype <- sctype_scores |> 
    left_join(type_short) |>
    group_by(type) |>
    arrange(-ncells) |> 
    mutate(type_rank = row_number(),
           cell_type_fine = paste0(cell_type_broad, 
                                   ".",
                                   width = str_pad(type_rank, 
                                                   nchar(as.character(max(type_rank))),
                                                   side = "left",
                                                   pad = "0")),
           cell_type_broad = factor(cell_type_broad, levels = cell_type_levels),
    )

## add levels to cell type fine
fine_levels <- sctype |> arrange(cell_type_broad) |> pull(cell_type_fine)
sctype <- sctype|>
    mutate(cell_type_fine = factor(cell_type_fine, levels = fine_levels)) |>
    arrange(cell_type_fine)

write.csv(sctype, file = here(data_dir, paste0("sctype_prelim.csv")), row.names = FALSE)
save(sctype, file = here(data_dir, paste0("sctype_prelim.Rdata")))

# slurmjobs::job_single('13_sctype_final', create_shell = TRUE, memory = '25G', command = "Rscript 13_sctype_final.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
