## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/b8719fbcdde7b14be267a396c1dd743cdc034649/code/09_HD_cell_level/05_banksy_embedding.R

library(here)
library(SpatialExperiment)
library(sessioninfo)
library(Banksy)
library(harmony)
library(cowplot)
library(scater)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'spe_banksy.rds'
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
plot_dir = here('plots', '09_HD_cell_level', 'no_secondary', 'banksy')

random_seed = 0
buffer_prop = 0.5
lambda = 0.2

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), showWarnings = FALSE)
set.seed(random_seed)

#   Load and subset to SVGs to avoid exceeding maximum number
#   of rows in a data table created internally. See
#   https://github.com/prabhakarlab/Banksy/issues/38#issuecomment-2310220881 and
#   https://github.com/prabhakarlab/Banksy_py/issues/12#issuecomment-2268114768
spe = readRDS(spe_path)
spe = spe[readLines(svg_path),]
spe$exclude_overlapping = FALSE

################################################################################
#   Stagger spatial coordinates to fit each sample in a unique range
################################################################################

message(Sys.time(), ' | Staggering spatial coordinates')
coords = spatialCoords(spe) |>
    as_tibble() |>
    mutate(sample_id = factor(spe$sample_id)) |>
    rename(sdimx = pxl_col_in_fullres, sdimy = pxl_row_in_fullres)

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

################################################################################
#   Compute Banksy embedding
################################################################################

message(Sys.time(), ' | Running computeBanksy on full dataset')
spe = computeBanksy(
    spe, assay_name = "logcounts", compute_agf = TRUE, seed = random_seed
)

message(Sys.time(), ' | Running PCA on embedding')
spe = runBanksyPCA(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

################################################################################
#   Run Harmony on embedding, with UMAP before and after
################################################################################

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

################################################################################
#   Explore effect of Harmony on UMAP
################################################################################

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
spe$batch_num = paste(
    'Batch',
    sample_info$batch_num[
        match(spe$sample_id, sample_info$tissue_id)
    ]
)

#   All samples together, colored by sample ID and batch number (separate plots)
for (color_var in c('sample_id', 'batch_num')) {
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
        point_alpha = 0.5, color_by = 'batch_num'
    ) +
        facet_wrap(~ spe$sample_id, nrow = 1) +
        theme_bw(base_size = 18) +
        guides(
            color = guide_legend(override.aes = list(size = 4, alpha = 1))
        ),
    plotReducedDim(
        spe, "UMAP_HARMONY", point_size = 0.6, point_alpha = 0.5,
        color_by = 'batch_num'
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

message(Sys.time(), ' | Saving full SPE object')
saveRDS(spe, file = out_path)

session_info()