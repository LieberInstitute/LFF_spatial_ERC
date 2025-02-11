## Louise Huuki-Myers, January 2025
## Annotate preliminary clusters with SC Type

## load libraries
library("dplyr")
library("purrr")
library("Seurat")
library("HDF5Array")
library("HGNChelper")
library("openxlsx")
library("here")
library("sessioninfo")

## source sc-type functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "08_sctype_prelim")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "08_sctype_prelim")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load sce object
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_harmony"))
sce

## load cluster outputs
load(here("processed-data", "04_snRNA-seq", "07_cluster_sn_prelim", "walktrap_snn_k%02d_clusters_prelim.Rdata"), verbose = TRUE)
# clusters

sce$snn_k10 <- sprintf("k10c%02d", clusters)

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

## TODO unnest
cluster_cols <- c("snn_k10")
names(cluster_cols) <- cluster_cols

sctype_scores <- map(cluster_cols, function(cl_col){
    
    clusters <- sort(unique(sce[[cl_col]]))
    
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
    
    write.csv(sctype_scores, file = here(data_dir, paste0("sctype_scores_", cl_col,".csv")))
    
    return(sctype_scores)
})

# slurmjobs::job_single('08_sctype_prelim', create_shell = TRUE, memory = '25G', command = "Rscript 08_sctype_prelim.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
