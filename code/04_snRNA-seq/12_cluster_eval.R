## Louise Huuki-Myers, February 2024
## Evaluate different clustering results with Silhouette width

library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("bluster")
library("cluster")
library("ggrepel")

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

message("number of clusters")
map_int(clusters, max)

message("number clusters n<100")
map_int(clusters, ~sum(table(.x) < 100))

#### Calculate Silhouette Score ####
message(Sys.time(), " - Calc Silhouette score")
sil.approx <- map(clusters, ~approxSilhouette(reducedDim(sce, "HARMONY"), clusters=.x))

message(Sys.time(), " - done...saving")
# save(sil.approx, file = here(data_dir, "sil.approx.Rdata"))

sil.approx[[1]]

## calc mean silhouette score
silhouette_mean <- map_dbl(sil.approx, ~mean(.x[, "width"]))

cluster_eval <- tibble(k = as.integer(gsub("k", "", colnames(clusters))),
                 n_clus = map_int(clusters, max),
                 silh_score = silhouette_mean)

cluster_eval |> arrange(-silh_score)

## k value with the highest silhouette score
optimal_k <-parse_number(names(silhouette_mean)[which.max(silhouette_mean)])
optimal_k
# [1] "k13"

cluster_eval <- cluster_eval |> mutate(silh_max = optimal_k == k)

write_csv(cluster_eval,  file = here(data_dir, "cluster_eval.csv"))


#### plot  ####
## k vs. n cluster
cluster_eval_k_vs_n <- ggplot(cluster_eval, aes(x = k, y = n_clus)) + ## add color by sil_score?
    geom_point() +
    geom_line() +
    geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
    theme_bw()

ggsave(cluster_eval_k_vs_n, filename = here(plot_dir, "cluster_eval_k_vs_n.png"))

## silhouette score vs. k
cluster_eval_k_vs_silh <- ggplot(cluster_eval, aes(x = k, y = silh_score)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
    theme_bw()  +
    labs(y = "Mean silhouette score")

ggsave(cluster_eval_k_vs_silh, filename = here(plot_dir, "cluster_eval_k_vs_silh.png"))

## silhouette score vs. n cluster
cluster_eval_n_vs_silh <- ggplot(cluster_eval, aes(x = n_clus, y = silh_score)) +
    geom_point(aes(color = k)) +
    geom_text_repel(aes(label = k)) +
    theme_bw() +
    labs(y = "Mean silhouette score")

ggsave(cluster_eval_n_vs_silh, filename = here(plot_dir, "cluster_eval_n_vs_silh.png"))

#### Silh distributions ####

sil.approx.all <- map2_dfr(sil.approx, names(sil.approx), ~as.data.frame(.x) |> mutate(k = .y))

sil_boxplot <- sil.approx.all |>
    ggplot(aes(x = k, y = width)) +
    geom_boxplot() +
    theme_bw()  +
    labs(y = "silhouette width")

ggsave(sil_boxplot, filename = here(plot_dir, "sil_boxplot.png"))


#### concordance matrix ####

k_combn <- combn(colnames(clusters), 2)

rand_scores <-  map2_dfr(k_combn[1,], k_combn[2,], ~list(kx = .x, ky = .y, rand = pairwiseRand(clusters[[.x]], clusters[[.y]], mode = "index")))

rand_scores |> arrange(rand)

rand_scores |>
    filter(kx == "k13" | ky == "k13") |>
    arrange(-rand)

rand_tile <- rand_scores |>
    ggplot(aes(kx, ky, fill = rand)) +
    geom_tile() +
    theme_bw() +
    viridis::scale_fill_viridis() +
    scale_x_discrete(position = "top") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggsave(rand_tile, filename = here(plot_dir, "rand_tile.png"))

# slurmjobs::job_single('12_cluster_eval', create_shell = TRUE, memory = '5G', command = "12_cluster_eval.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

