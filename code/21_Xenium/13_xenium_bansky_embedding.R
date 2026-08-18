## Louise Huuki-Myers, May 2026
## Run Bansky on ERC Xenium data to get SpDs
## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/b8719fbcdde7b14be267a396c1dd743cdc034649/code/09_HD_cell_level/05_banksy_embedding.R

#### Set Up ####

library("here")
library("SpatialExperiment")
library("sessioninfo")
library("Banksy")
library("harmony")
library("cowplot")
library("scater")
library("tidyverse")
library("qs2")

library("spatialLIBD")
library("scDotPlot")

data_dir <- here("processed-data", "21_Xenium", "13_xenium_bansky_embedding")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "13_xenium_bansky_embedding")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data #### 
message(Sys.time(), " - Load spe data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message("n cells: ", ncol(spe))

#### set up params ####
random_seed = 514
buffer_prop = 0.5
lambda = 0.8 ## 0.8 finds spatial domain vs. 0.2 cell type resolution
k_geom <- c(25, 50)

set.seed(random_seed)


####   Stagger spatial coordinates to fit each sample in a unique range ####

orginal_coords <- spatialCoords(spe)

message(Sys.time(), ' | Staggering spatial coordinates')
coords = spatialCoords(spe) |>
    as_tibble() |>
    mutate(sample_id = factor(spe$sample_id)) |>
    rename(sdimx = x_centroid, sdimy = y_centroid)

#   Find a range of X values slightly larger than any particular sample
x_size = coords |>
    group_by(sample_id) |>
    summarize(x_diff = max(sdimx) - min(sdimx)) |>
    pull(x_diff) |>
    max()
x_size = (1 + buffer_prop) * x_size

#   Separate samples by placing each sample into the same Y range and adjacent
#   X ranges
spatialCoords(spe) = coords |>
    group_by(sample_id) |>
    mutate(
        sdimx = sdimx - min(sdimx) + x_size * (match(cur_group()$sample_id, unique(spe$sample_id)) - 1),
        sdimy = sdimy - min(sdimy)
    ) |>
    ungroup() |>
    select(sdimx, sdimy) |>
    as.matrix()

#### Compute Banksy embedding ####

message(Sys.time(), ' | Running computeBanksy on full dataset')
spe = computeBanksy(
    spe, assay_name = "logcounts", compute_agf = TRUE, seed = random_seed, k_geom = k_geom
)

message(Sys.time(), ' | Running PCA on embedding')
spe = runBanksyPCA(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

####   Run Harmony on embedding, with UMAP before and after ####

message(Sys.time(), ' | Running UMAP on embedding')
spe = runBanksyUMAP(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

message(Sys.time(), " | Running Harmony...")
reducedDims(spe)$PCA = reducedDims(spe)[[sprintf('PCA_M1_lam%s', lambda)]]
reducedDims(spe)[[sprintf('PCA_M1_lam%s', lambda)]] = NULL
pdf(file.path(plot_dir, "harmony_convergence.pdf"))
spe = RunHarmony(
    spe, group.by.vars = "sample_id", plot_convergence = TRUE
)
dev.off()

message(Sys.time(), ' | Running UMAP on Harmony-corrected embedding')
spe = runBanksyUMAP(
    spe,  dimred = "HARMONY", use_agf = TRUE, lambda = lambda,
    seed = random_seed
)

####   Explore effect of Harmony on UMAP ####

#   All samples together, colored by sample ID and batch number (separate plots)
for (color_var in c('sample_id', 'chip')) {
    p = plot_grid(
        plotReducedDim(
            spe, sprintf("UMAP_M1_lam%s", lambda), point_size = 0.6,
            point_alpha = 0.5, color_by = color_var
        ) +
            theme_classic(base_size = 18) +
            theme(legend.position = "none"),
        plotReducedDim(
            spe, "UMAP_HARMONY", point_size = 0.6, point_alpha = 0.5,
            color_by = color_var
        ) +
            theme_classic(base_size = 18) +
            guides(
                color = guide_legend(override.aes = list(size = 4, alpha = 1))
            ),
        nrow = 1,
        rel_widths = c(1, 1.2)
    )
    png(
        file.path(
            plot_dir, sprintf('harmony_umap_%s_together.png', color_var)
        ),
        width = 1200, height = 600
    )
    print(p)
    dev.off()
}

#   Faceted by sample, colored by batch number
p = plot_grid(
    plotReducedDim(
        spe, sprintf("UMAP_M1_lam%s", lambda), point_size = 0.6,
        point_alpha = 0.5, color_by = 'chip'
    ) +
        facet_wrap(~ spe$sample_id, nrow = 1) +
        theme_bw(base_size = 18) +
        guides(
            color = guide_legend(override.aes = list(size = 4, alpha = 1))
        ),
    plotReducedDim(
        spe, "UMAP_HARMONY", point_size = 0.6, point_alpha = 0.5,
        color_by = 'chip'
    ) +
        facet_wrap(~ spe$sample_id, nrow = 1) +
        theme_bw(base_size = 18) +
        guides(
            color = guide_legend(override.aes = list(size = 4, alpha = 1))
        ),
    nrow = 2
)
png(
    file.path(plot_dir, 'harmony_umap_apart.png'),
    width = 1500, height = 750
)
print(p)
dev.off()

#### Bansky clustering ####
# swap back to OG coords
spatialCoords(spe) <- orginal_coords

cluster_colors <- c("#f62062",
                    "#f45e28",
                    "#cf9800",
                    "#608d00",
                    "#01e090",
                    "#3ed9e6",
                    "#0064ca",
                    "#9215a3",
                    "#ff8ee2",
                    "black",
                    "grey",
                    "brown")

## app data check (k=9 specific spot check, kept for reference)
# spe_app <- qs_read(here("code", "23_spatialLIBD_app_Xenium", "spe_xenium_app.qs2"))
# table(spe_app[, colnames(spe)]$clust_HARMONY_kmeans9)
# table(spe_app[, colnames(spe)]$clust_HARMONY_kmeans9, spe$clust_HARMONY_kmeans9)
# colData(spe)[,duplicated(colnames(colData(spe)))] <- NULL

for (k in 9:12) {
    message(Sys.time(), sprintf(' | Performing clustering k=%d', k))
    spe = clusterBanksy(
        spe, 
        use_agf = TRUE, 
        lambda = lambda,
        seed = random_seed,
        algo = "kmeans", 
        # resolution = res, 
        dimred = "HARMONY",
        kmeans.centers = k
    )

    clust_var <- sprintf("clust_HARMONY_kmeans%d", k)
    print(table(colData(spe)[[clust_var]]))

    map(unique(spe$BrNum), function(samp){
        vis_clus_class <- spatialLIBD::vis_clus(spe,
                                                sampleid = samp,
                                                clustervar = clust_var,
                                                datatype = "Xenium",
                                                point_size = 1.5,
                                                # alpha = 0.5,
                                                colors = cluster_colors,
                                                guide_point_size = 3)

        ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_bansky_k%d_%s.png", k, samp)), width = 12)

    })
}

table(spe$clust_HARMONY_kmeans9, spe$clust_HARMONY_kmeans12)

#### Save preliminary SPE with Bansky embeddings + clusters ####
message(Sys.time(), " - Saving preliminary SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_bansky_prelim.qs2"))

# spe <- qs_read(here(data_dir, "spe_xenium_bansky_prelim.qs2"))


# slurmjobs::job_single('13_xenium_bansky_embedding', create_shell = TRUE, memory = '100G', command = "Rscript 13_xenium_bansky_embedding.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
