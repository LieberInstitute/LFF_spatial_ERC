## Louise Huuki-Myers, February 2024
## Evaluate different clustering results with Silhouette width

library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("bluster")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "12_cluster_eval")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "12_cluster_eval")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), "- load reprocessed sce")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_reprocess"))
sce

#### load cluster outputs ####
message(Sys.time(), " - Load cluster data")
cluster_fn <- list.files(here("processed-data", "04_snRNA-seq", "11_cluster_sn"), full.names = TRUE)  
names(cluster_fn) <- ss(basename(cluster_fn), "_", 3)

clusters <- map_dfc(cluster_fn, ~get(load(.x)))
head(clusters)
dim(clusters)

map_int(clusters, max)
# k10 k15 k20 
# 62  28  33 

map_int(clusters, ~sum(table(.x) < 100))
# k10 k15 k20 
# 12   4   1 


#### Calculate Silhouette Score ####
message(Sys.time(), " - Calc Distantce matrix")
dist_matrix <- dist(reducedDim(sce, "HARMONY"))

message(Sys.time(), " - Calc Silhouette score")
silhouette_score <- map(clusters, ~silhouette(.x, dist_matrix))
save(silhouette_score, file = "silhouette_score.Rdata")

## calc mean silhouette score
silhouette_mean <- map(silhouette_score, ~mean(.x[, "sil_width"]))

cluster_eval <- tibble(k = as.integer(gsub("k", "", colnames(clusters))),
                 n_clus = map_int(clusters, max),
                 silh_score = silhouette_mean)

write_csv(cluster_eval,  file = here(data_dir, "cluster_eval.csv"))

## k value with the highest silhouette score
optimal_k <- names(silhouette_mean)[which.max(silhouette_mean)]
optimal_k
# [1] 20


#### plot  ####
## k vs. n cluster
cluster_eval_k_vs_n <- ggplot(cluster_eval, aes(x = k, y = n_clus)) + ## add color by sil_score?
    geom_point() +
    geom_line() +
    theme_bw()

ggsave(cluster_eval_k_vs_n, filename = here(plot_dir, "cluster_eval_k_vs_n.png"))

## silhouette score vs. k
cluster_eval_k_vs_silh <- ggplot(cluster_eval, aes(x = k, y = silh_score)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave(cluster_eval_k_vs_silh, filename = here(plot_dir, "cluster_eval_k_vs_silh.png"))

## silhouette score vs. n cluster
cluster_eval_n_vs_silh <- ggplot(cluster_eval, aes(x = n_clus, y = silh_score)) +
    geom_point() +
    # geom_line() +
    theme_bw()

ggsave(cluster_eval_n_vs_silh, filename = here(plot_dir, "cluster_eval_n_vs_silh.png"))


# slurmjobs::job_single('12_cluster_eval', create_shell = TRUE, memory = '5G', command = "12_cluster_eval.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

