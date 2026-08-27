## Louise Huuki-Myers, May 2026
## Annotate Bansky spatial domains (SpX) for ERC Xenium data
## Continues from 13_xenium_bansky_embedding.R (reads in spe_xenium_bansky_prelim.qs2)

#### Set Up ####

library("here")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("qs2")

library("spatialLIBD")
library("scDotPlot")
library("ComplexHeatmap")

data_dir <- here("processed-data", "21_Xenium", "13_xenium_bansky_embedding")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "13_xenium_bansky_embedding")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load prelim data ####
message(Sys.time(), " - Load prelim spe data")
spe <- qs_read(here(data_dir, "spe_xenium_bansky_prelim.qs2"))

message("n cells: ", ncol(spe))

#### Annotate bansky clus ####

load(here("processed-data","SpX_colors.Rdata"), verbose = TRUE)

#### Annotate bansky clus (k11) ####
## banksy_clustering_annotation_k11_summary.xlsx was built against the
## already-relabeled clust_HARMONY_kmeans11 values ("xSp11D01" etc.), so this
## one matches directly, unlike the k9/SpX crosswalk above
xSpD_anno_table <- readxl::read_xlsx(here(
    "processed-data", "21_Xenium", "13.1_xenium_compare_bansky_clusters", "banksy_clustering_annotation_k11_summary.xlsx"
)) |>
    mutate(xSpD = fct_reorder(paste0("x", anno), order),
           xSpD_k11 = fct_reorder(paste0(anno, "~", cluster), order))

spe$xSpD <- xSpD_anno_table$xSpD[match(spe$clust_HARMONY_kmeans11, xSpD_anno_table$cluster)]
spe$xSpD_k11 <- xSpD_anno_table$xSpD[match(spe$clust_HARMONY_kmeans11, xSpD_anno_table$cluster)]

table(spe$xSpD, spe$clust_HARMONY_kmeans11)
stopifnot(!anyNA(spe$xSpD))

## plain domain name (no cluster id), used for plots/heatmaps/color-matching -
## same role vSpD_anno plays alongside vSpD in the Visium scripts
spe$xSpD_anno <- xSpD_anno_table$anno[match(spe$clust_HARMONY_kmeans11, xSpD_anno_table$cluster)]
stopifnot(!anyNA(spe$xSpD_anno))

cluster_data <- as.data.frame(colData(spe)[,c("sample_id", "Barcode",
                                              "clust_HARMONY_kmeans9",
                                              "clust_HARMONY_kmeans10",
                                              "clust_HARMONY_kmeans11",
                                              "clust_HARMONY_kmeans12",
                                              "xSpD",
                                              "xSpD_k11")])
write.csv(cluster_data, file = here(data_dir, "Xenium_bansky_cluster_data.csv"))

#### Cell type composition (RCTD label transfer) ####
## Loaded early because cell-type composition feeds into the annotation
message(Sys.time(), " - Load RCTD data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD", "rctd_results_xenium.qs2"))

spe$cell_id <- colnames(spe)
colData(spe) <- cbind(colData(spe), rctd_data@results$results_df[colnames(spe), ])

spe$cell_type_anno <- spe$first_type
spe$cell_type_broad <- factor(
    gsub("\\..*?$", "", spe$cell_type_anno),
    levels = c("Astro", "Macro", "Micro", "Oligo", "OPC", "Vasc", "Excit", "Inhib")
)

table(spe$cell_type_anno)
table(spe$cell_type_broad)



#### Annotate Regions ####

region_lookup <- c(
    "Vasc"     = "Vasc",
    "Vasc.im"  = "Vasc",
    "Vasc.adj" = "Vasc",
    "GL"       = "GM_L1",
    "L2"       = "GM_L2_6",
    "Inhib"    = "GM_L2_6",
    "LD"       = "GM_L2_6",
    "L5"       = "GM_L2_6",
    "L6"       = "GM_L2_6",
    "WMuf"     = "WM",
    "WMd"      = "WM"
)

spe$region <- factor(region_lookup[spe$xSpD_anno], levels = c("Vasc", "GM_L1", "GM_L2_6", "WM"))

table(spe$xSpD, spe$region)

region_colors = c('Vasc' = "#E05AD2",
                  'GM_L1' = "#67a75f",
                  'GM_L2_6' = "#3446BD",
                  'WM'= "gray")

metadata(spe)$region_colors <- region_colors

#### Save SPE with bansky data ####
message(Sys.time(), " - Saving SPE object")

metadata(spe)$SpX_colors <- SpX_colors
metadata(spe)$SpX_colors <- SpX_colors

qs2::qs_save(spe, here(data_dir, "spe_xenium_bansky.qs2"))

# spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

#### Plot Vis clus xSpD ####

map(unique(spe$BrNum), function(samp){
    vis_clus_class <- spatialLIBD::vis_clus(spe,
                                            sampleid = samp,
                                            clustervar = "xSpD",
                                            datatype = "Xenium",
                                            point_size = 1.5,
                                            # alpha = 0.5,
                                            colors = SpX_colors,
                                            guide_point_size = 3)
    
    ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_SpX_%s.png", samp)), width = 12)
    
})

map(unique(spe$BrNum), function(samp){
    vis_clus_class <- spatialLIBD::vis_clus(spe,
                                            sampleid = samp,
                                            clustervar = "region",
                                            datatype = "Xenium",
                                            point_size = 1.5,
                                            # alpha = 0.5,
                                            colors = region_colors,
                                            guide_point_size = 3)
    
    ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_region_%s.png", samp)), width = 12)
    
})


source(here("code", "utils", "my_plot_reduced_dim.R"))

## previously plotted clsuters, could flatten loop
walk2(c("xSpD"), list(SpX_colors),  ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP_HARMONY", my_var = .x, var_type = "cat", color_pal = .y))


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

map(spd_marker_list, ~all(.x %in% rownames(spe)))

DeconvoBuddies::plot_marker_express_List(sce = spe, 
                                         spd_marker_list, 
                                         color_pal = SpX_colors,
                                         cellType_col = "xSpD",
                                         pdf = here(plot_dir, "xenium_bansky_xSpD_SpD_markers.pdf"),
                                         gene_name_col = "gene_name")

## APOE expression
spx_APOE <- DeconvoBuddies::plot_gene_express(sce = spe, 
                                  genes = "APOE",
                                  category = "xSpD",
                                  color_pal = SpX_colors
                                  )

ggsave(spx_APOE, filename = here(plot_dir, "xenium_bansky_xSpD_APOE.png"))

apoe_mean <- as.data.frame(map_dbl(rafalib::splitit(spe$xSpD_anno), ~mean(logcounts(spe)["APOE",.x])))
colnames(apoe_mean) <- "APOE_mean"

write.csv(apoe_mean, file = here(data_dir, "Xenium_xSpD_APOE_mean_logcount.csv"))

#### scDot plots ####

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

rowData(spe)$SpD_marker <- NULL
rowData(spe)$SpD_marker <- marker_stats_top$SpD[match(rownames(spe), marker_stats_top$gene)] 
table(rowData(spe)$SpD_marker)


pdf(here(plot_dir, "Xenium_bansky_xSpD_dotplot_SpD_markers.pdf"))
spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "SpD_marker",
              scale = TRUE,
              annoColors = list(xSpD_anno = SpX_colors,
                                SpD_marker = SpD_colors),
              clusterRows = FALSE,
              clusterColumns = TRUE,
              groupLegends = FALSE)

spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "SpD_marker",
              scale = TRUE,
              annoColors = list(xSpD_anno = SpX_colors,
                                SpD_marker = SpD_colors),
              clusterRows = FALSE,
              clusterColumns = FALSE,
              groupLegends = FALSE)
dev.off()

#### Cell Type vs. xSpD heatmap ####

load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)

cell_v_xSpD <- table(spe[, spe$spot_class == "singlet"]$xSpD, spe[, spe$spot_class == "singlet"]$cell_type_anno)

## proportion cell type (xSpD rows sum to 1)
cell_prop_v_xSpD <- sweep(cell_v_xSpD, 1, rowSums(cell_v_xSpD), FUN = "/")
rowSums(cell_prop_v_xSpD)

## proportion xSpD (cell type cols sum to 1)
cell_v_xSpD_prop <- sweep(cell_v_xSpD, 2, colSums(cell_v_xSpD), FUN = "/")
colSums(cell_v_xSpD_prop)

## save & order proportion data 

cell_v_xSpD_prop_long <- cell_v_xSpD |> 
    reshape2::melt() |>
    rename(xSpD = Var1, cell_type_anno = Var2, n_cell = value) |>
    left_join(cell_prop_v_xSpD |> 
                  reshape2::melt() |>
                  rename(xSpD = Var1, cell_type_anno = Var2, prop_cell_type = value)) |>
    left_join(cell_v_xSpD_prop |> 
                  reshape2::melt() |>
                  rename(xSpD = Var1, cell_type_anno = Var2, prop_xSpD = value)) |>
    as_tibble()

write.csv(cell_v_xSpD_prop_long, file = here(data_dir, "cell_v_xSpD_prop_long.csv"))

cell_v_xSpD_max <- cell_v_xSpD_prop_long |> 
    group_by(cell_type_anno) |>
    slice_max(n_cell) |> 
    arrange(xSpD, cell_type_anno) |>
    print(n= 38)


cell_v_xSpD_prop_long |> 
    filter(cell_type_anno == "Oligo.3") |>
    arrange(-n_cell)

## reorder cell types by max xSpD
cell_v_xSpD <- cell_v_xSpD[, cell_v_xSpD_max$cell_type_anno]
cell_prop_v_xSpD <- cell_prop_v_xSpD[, cell_v_xSpD_max$cell_type_anno]
cell_v_xSpD_prop <- cell_v_xSpD_prop[, cell_v_xSpD_max$cell_type_anno]

## create annotations 
xSpD_row_ha <- rowAnnotation(
    xSpD = rownames(cell_v_xSpD),
    col = list(xSpD = SpX_colors),
    show_legend = FALSE
)

all(colnames(cell_v_xSpD) %in% names(cell_type_colors$anno))

cell_type_col_ha <- HeatmapAnnotation(
    cell_type = colnames(cell_v_xSpD),
    col = list(cell_type = cell_type_colors$anno),
    show_legend = FALSE
)


## PLOT HEATMAPS
pdf(here(plot_dir, "Xenium_bansky_xSpD_v_cell_type_heatmap.pdf"), width = 10)

## counts
ComplexHeatmap::Heatmap(cell_v_xSpD,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_v_xSpD,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)

## proportion cell type (xSpD row sum to 1)
ComplexHeatmap::Heatmap(cell_prop_v_xSpD,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,  ## cluster xSpD
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_prop_v_xSpD,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = TRUE, ## Clsuter cell types
                        cluster_rows = TRUE,  ## cluster xSpD
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_prop_v_xSpD,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)

## proportion xSpD (cell type cols sum to 1)
ComplexHeatmap::Heatmap(cell_v_xSpD_prop,
                        name = "prop xSpD\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = TRUE,
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_v_xSpD_prop,
                        name = "prop xSpD\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha, 
                        bottom_annotation = cell_type_col_ha)
dev.off()


## Add cell type hierarchical clusters 

# ## Load Hierarchical clustering
# load(here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster","sn_subcluster_hierarchical_cluster.Rdata"), verbose = TRUE)

## PLOT HEATMAPS
pdf(here(plot_dir, "Xenium_bansky_xSpD_v_cell_type_heatmap_broad.pdf"), width = 10)

## proportion xSpD (cell type cols sum to 1)
ComplexHeatmap::Heatmap(cell_v_xSpD_prop,
                        name = "prop xSpD\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        # cluster_columns = FALSE,
                        column_split = gsub("\\..*", "", colnames(cell_v_xSpD_prop)),
                        # cluster_columns = cluster_within_group(cell_v_xSpD_prop, gsub("\\..*", "", colnames(cell_v_xSpD_prop))),
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha,
                        bottom_annotation = cell_type_col_ha
                        )

ComplexHeatmap::Heatmap(cell_v_xSpD,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        column_split = gsub("\\..*", "", colnames(cell_v_xSpD)),
                        cluster_rows = FALSE, 
                        left_annotation = xSpD_row_ha,
                        bottom_annotation = cell_type_col_ha
                        )
dev.off()


#### Cell Type vs. region heatmap ####

cell_v_region <- table(spe[, spe$spot_class == "singlet"]$region, spe[, spe$spot_class == "singlet"]$cell_type_anno)

## proportion cell type (region rows sum to 1)
cell_prop_v_region <- sweep(cell_v_region, 1, rowSums(cell_v_region), FUN = "/")
rowSums(cell_prop_v_region)

## proportion region (cell type cols sum to 1)
cell_v_region_prop <- sweep(cell_v_region, 2, colSums(cell_v_region), FUN = "/")
colSums(cell_v_region_prop)

## save & order proportion data 

cell_v_region_prop_long <- cell_v_region |> 
    reshape2::melt() |>
    rename(region = Var1, cell_type_anno = Var2, n_cell = value) |>
    left_join(cell_prop_v_region |> 
                  reshape2::melt() |>
                  rename(region = Var1, cell_type_anno = Var2, prop_cell_type = value)) |>
    left_join(cell_v_region_prop |> 
                  reshape2::melt() |>
                  rename(region = Var1, cell_type_anno = Var2, prop_region = value)) |>
    as_tibble()

write.csv(cell_v_region_prop_long, file = here(data_dir, "cell_v_region_prop_long.csv"))

cell_v_region_max <- cell_v_region_prop_long |> 
    group_by(cell_type_anno) |>
    slice_max(n_cell) |> 
    arrange(region, cell_type_anno) 

cell_v_region_max |>
    print(n= 38)


cell_v_region_prop_long |> 
    filter(cell_type_anno == "Oligo.3") |>
    arrange(-n_cell)

## reorder cell types by max region
cell_v_region <- cell_v_region[, cell_v_region_max$cell_type_anno]
cell_prop_v_region <- cell_prop_v_region[, cell_v_region_max$cell_type_anno]
cell_v_region_prop <- cell_v_region_prop[, cell_v_region_max$cell_type_anno]

## create annotations 
region_row_ha <- rowAnnotation(
    region = rownames(cell_v_region),
    col = list(region = region_colors),
    show_legend = FALSE
)

all(colnames(cell_v_region) %in% names(cell_type_colors$anno))

cell_type_col_ha <- HeatmapAnnotation(
    cell_type = colnames(cell_v_region),
    col = list(cell_type = cell_type_colors$anno),
    show_legend = FALSE
)


## PLOT HEATMAPS
pdf(here(plot_dir, "Xenium_bansky_region_v_cell_type_heatmap.pdf"), width = 10)

## counts
ComplexHeatmap::Heatmap(cell_v_region,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_rows = FALSE, 
                        left_annotation = region_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_v_region,
                        name = "n singlet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = region_row_ha, 
                        bottom_annotation = cell_type_col_ha)

## proportion cell type (region row sum to 1)
ComplexHeatmap::Heatmap(cell_prop_v_region,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,  ## cluster region
                        left_annotation = region_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_prop_v_region,
                        name = "prop cell type\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = region_row_ha, 
                        bottom_annotation = cell_type_col_ha)

## proportion region (cell type cols sum to 1)
ComplexHeatmap::Heatmap(cell_v_region_prop,
                        name = "prop region\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = TRUE,
                        cluster_rows = FALSE, 
                        left_annotation = region_row_ha, 
                        bottom_annotation = cell_type_col_ha)

ComplexHeatmap::Heatmap(cell_v_region_prop,
                        name = "prop region\nsinglet cells",
                        col = c("black", viridisLite::plasma(100)),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE, 
                        left_annotation = region_row_ha, 
                        bottom_annotation = cell_type_col_ha)
dev.off()



# slurmjobs::job_single('13.5_xenium_bansky_annotate', create_shell = TRUE, memory = '50G', command = "Rscript 13.5_xenium_bansky_annotate.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()