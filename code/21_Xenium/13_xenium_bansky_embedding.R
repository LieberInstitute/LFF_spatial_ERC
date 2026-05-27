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

message(Sys.time(), ' | Performing clustering k=9')
spe = clusterBanksy(
    spe, 
    use_agf = TRUE, 
    lambda = lambda,
    seed = random_seed,
    algo = "kmeans", 
    # resolution = res, 
    dimred = "HARMONY",
    kmeans.centers = 9
)

# spe_app <- qs_read(here("code", "23_spatialLIBD_app_Xenium", "spe_xenium_app.qs2"))
# 
# table(spe_app[, colnames(spe)]$clust_HARMONY_kmeans9)
# table(spe_app[, colnames(spe)]$clust_HARMONY_kmeans9, spe$clust_HARMONY_kmeans9)

spe$clust_HARMONY_kmeans9 <- spe_app[, colnames(spe)]$clust_HARMONY_kmeans9
spe$clust_HARMONY_kmeans12 <- spe_app[, colnames(spe)]$clust_HARMONY_kmeans12

table(spe$clust_HARMONY_kmeans9)

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

# samp <- "Br1556"

map(unique(spe$BrNum), function(samp){
    vis_clus_class <- spatialLIBD::vis_clus(spe,
                                            sampleid = samp,
                                            clustervar = "clust_HARMONY_kmeans9",
                                            datatype = "Xenium",
                                            point_size = 1.5,
                                            # alpha = 0.5,
                                            colors = cluster_colors,
                                            guide_point_size = 3)
    
    ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_bansky_k9_%s.png", samp)), width = 12)
    
})

message(Sys.time(), ' | Performing clustering k=12')
spe = clusterBanksy(
    spe, 
    use_agf = TRUE, 
    lambda = lambda,
    seed = random_seed,
    algo = "kmeans", 
    # resolution = res, 
    dimred = "HARMONY",
    kmeans.centers = 12
)

map(unique(spe$BrNum), function(samp){
    vis_clus_class <- spatialLIBD::vis_clus(spe,
                                            sampleid = samp,
                                            clustervar = "clust_HARMONY_kmeans12",
                                            datatype = "Xenium",
                                            point_size = 1.5,
                                            # alpha = 0.5,
                                            colors = cluster_colors,
                                            guide_point_size = 3)
    
    ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_bansky_k12_%s.png", samp)), width = 12)
    
})

table(spe$clust_HARMONY_kmeans12)
table(spe$clust_HARMONY_kmeans9, spe$clust_HARMONY_kmeans12)

colData(spe)

#### Annotate bansky clus ####

SpX_colors = c('Vasc~SpX3' = "#E05AD2",
               'L1~SpX6' = "#9AA7FE",
               'L1~SpX7' = "#0220DE",
               'L2.3~SpX4' = "#FEAF16",
               'Inhib~SpX5' = "#C82100",
               'L5~SpX1' = "#16FF32",
               'L6~SpX9' = "#178C6D",
               'WMtz~SpX8' = "grey",
               'WM~SpX2'= "#581009")

cluster_anno <- readxl::read_xlsx(here("processed-data", "21_Xenium", "Bansky_cluster_notes.xlsx")) |>
    filter(!grepl("uf", Visium)) |>
    mutate(SpX = factor(paste0(Visium , "~SpX", Xenium_k9), levels = names(SpX_colors)))

# load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
# SpD_colors
# Vasc~Sp09D08    L1~Sp09D05  L2.3~Sp09D01    LD~Sp09D02 Inhib~Sp09D09    L5~Sp09D03    L6~Sp09D04 WM.uf~Sp09D07
# "#E05AD2"     "#0220DE"     "#FEAF16"     "#00BCF9"     "#C82100"     "#16FF32"     "#178C6D"     "#E4E1E3"
# WM~Sp09D06
# "#581009"

spe$SpX <- cluster_anno$SpX[match(spe$clust_HARMONY_kmeans9, cluster_anno$Xenium_k9)]

table(spe$SpX, spe$clust_HARMONY_kmeans9)

#### Add RCTD data ###

message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

rctd_data@results$results_df <- rctd_data@results$results_df[colnames(spe),]

identical(colnames(spe), rownames(rctd_data@results$results_df))

spe$cell_id <- colnames(spe)
colData(spe) <- cbind(colData(spe), rctd_data@results$results_df)

colData(spe) <- cbind(colData(spe), rctd_data@results$results_df)

spe$cell_type_anno <- spe$first_type
spe$cell_type_broad <- factor(gsub("\\..*?$", "", spe$cell_type_anno), levels = c("Astro", "Macro", "Micro","Oligo", "OPC","Vasc","Excit","Inhib"))
table(spe$cell_type_broad)

table(spe$SpX, spe$cell_type_broad)
table(spe$SpX, spe$spot_class)

#### Save SPE with bansky data ####
message(Sys.time(), " - Saving SPE object")

metadata(spe)$SpX_colors <- SpX_colors

qs2::qs_save(spe, here(data_dir, "spe_xenium_bansky.qs2"))

# spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

#### Plot Vis clus SpX ####

# pdf(here(plot_dir, "Bansky_k9_cluster_v_cell_type_broad.pdf"))
# ComplexHeatmap::Heatmap(table(spe$clust_HARMONY_kmeans9, spe$cell_type_broad))
# dev.off()

map(unique(spe$BrNum), function(samp){
    vis_clus_class <- spatialLIBD::vis_clus(spe,
                                            sampleid = samp,
                                            clustervar = "SpX",
                                            datatype = "Xenium",
                                            point_size = 1.5,
                                            # alpha = 0.5,
                                            colors = SpX_colors,
                                            guide_point_size = 3)
    
    ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_SpX_%s.png", samp)), width = 12)
    
})


source(here("code", "utils", "my_plot_reduced_dim.R"))


walk2(c("clust_HARMONY_kmeans9", "SpX"), list(cluster_colors, SpX_colors),  ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP_HARMONY", my_var = .x, var_type = "cat", color_pal = .y))


#### SpD gene expression ####

load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_SpD.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(gene %in% rownames(spe), MeanRatio >1) |>
    group_by(cellType.target) |>
    arrange(MeanRatio.rank) |>
    slice(1:5) |>
    select(gene, MeanRatio.rank, SpD = cellType.target) |>
    arrange(SpD)

spd_marker_list <- map(rafalib::splitit(marker_stats_top$SpD), ~marker_stats_top$gene[.x])
spd_marker_list <- spd_marker_list[map_int(spd_marker_list, length) > 0]

map(lit_markers_list, ~all(.x %in% rownames(spe)))
map(lit_markers_list, ~all(.x %in% rownames(spe)))

DeconvoBuddies::plot_marker_express_List(sce = spe, 
                                         spd_marker_list, 
                                         color_pal = cluster_colors,
                                         cellType_col = "clust_HARMONY_kmeans9",
                                         pdf = here(plot_dir, "xenium_bansky_cluster_k9_SpD_markers.pdf"),
                                         gene_name_col = "Symbol")

DeconvoBuddies::plot_marker_express_List(sce = spe, 
                                         spd_marker_list, 
                                         color_pal = SpX_colors,
                                         cellType_col = "SpX",
                                         pdf = here(plot_dir, "xenium_bansky_SpX_SpD_markers.pdf"),
                                         gene_name_col = "Symbol")

#### scDot plots ####

rowData(spe)$SpD_marker <- NULL
rowData(spe)$SpD_marker <- marker_stats_top$SpD[match(rownames(spe), marker_stats_top$gene)] 
table(rowData(spe)$SpD_marker)

pdf(here(plot_dir, "Xenium_bansky_k9_dotplot_SpD_markers.pdf"))
spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "clust_HARMONY_kmeans9",
              groupAnno = "clust_HARMONY_kmeans9",
              featureAnno = "SpD_marker",
              scale = TRUE,
              annoColors = list(clust_HARMONY_kmeans9 = cluster_colors,
                                SpD_marker = SpD_colors),
              clusterRows = FALSE,
              clusterColumns = TRUE,
              groupLegends = FALSE)
dev.off()

pdf(here(plot_dir, "Xenium_bansky_SpX_dotplot_SpD_markers.pdf"))
spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "SpX",
              groupAnno = "SpX",
              featureAnno = "SpD_marker",
              scale = TRUE,
              annoColors = list(SpX = SpX_colors,
                                SpD_marker = SpD_colors),
              clusterRows = FALSE,
              clusterColumns = TRUE,
              groupLegends = FALSE)
dev.off()

#### Cell Type vs. SpX heatmap ####
library(ComplexHeatmap)

load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)

cell_v_SpX <- table(spe[, spe$spot_class == "singlet"]$SpX, spe[, spe$spot_class == "singlet"]$cell_type_anno)

## proportion cell type (SpX rows sum to 1)
cell_prop_v_SpX <- sweep(cell_v_SpX, 1, rowSums(cell_v_SpX), FUN = "/")
rowSums(cell_prop_v_SpX)

## proportion SpX (cell type cols sum to 1)
cell_v_SpX_prop <- sweep(cell_v_SpX, 2, colSums(cell_v_SpX), FUN = "/")
colSums(cell_v_SpX_prop)


SpX_row_ha <- rowAnnotation(
    SpX = rownames(cell_v_SpX),
    col = list(SpX = SpX_colors),
    show_legend = FALSE
)

all(colnames(cell_v_SpX) %in% names(cell_type_colors$anno))

cell_type_col_ha <- HeatmapAnnotation(
    cell_type = colnames(cell_v_SpX),
    col = list(cell_type = cell_type_colors$anno),
    show_legend = FALSE
)

## PLOT HEATMAPS
pdf(here(plot_dir, "Xenium_bansky_SpX_v_cell_type_heatmap.pdf"), width = 10)

## counts
ComplexHeatmap::Heatmap(cell_v_SpX,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_rows = FALSE, 
                        left_annotation = SpX_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_v_SpX,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = SpX_row_ha, 
                        bottom_annotation = cell_type_col_ha)

## proportion cell type (SpX row sum to 1)
ComplexHeatmap::Heatmap(cell_prop_v_SpX,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = TRUE, 
                        left_annotation = SpX_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_prop_v_SpX,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = SpX_row_ha, 
                        bottom_annotation = cell_type_col_ha)

## proportion SpX (cell type cols sum to 1)
ComplexHeatmap::Heatmap(cell_v_SpX_prop,
                        name = "prop SpX\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = TRUE,
                        cluster_rows = FALSE, 
                        left_annotation = SpX_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_v_SpX_prop,
                        name = "prop SpX\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = SpX_row_ha, 
                        bottom_annotation = cell_type_col_ha)
dev.off()


# slurmjobs::job_single('13_xenium_bansky_embedding', create_shell = TRUE, memory = '100G', command = "Rscript 13_xenium_bansky_embedding.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()