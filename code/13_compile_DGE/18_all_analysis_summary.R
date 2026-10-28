## Louise Huuki-Myers July 2026
## All analysis summary
##
## REWORK: t-stat heatmaps across the DE resolution cascade (Section 2.4-2.6
## of the preprint: spatial domain -> broad cell type -> fine subcluster),
## using logFC_Heatmap() with fill_stat = "vlmf_t" for a fixed anchor gene
## panel that recurs across all three levels. Xenium validation is shown as a
## 4th block: Xenium_cell_type_anno and Xenium_cell_type_anno_SpX (clusters =
## cell_type_anno x xSpD, stored as cluster_xSpD). Xenium_Oligo.3_Astro is
## loaded but has no heatmap yet.


#### set up ####
library("here")
library("tidyverse")
library("qs2")
library("sessioninfo")
library("patchwork")
library("ComplexHeatmap")
library("circlize")

data_dir <- here("processed-data", "13_compile_DGE", "18_all_analysis_summary")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "18_all_analysis_summary")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
# APOE_carrier_colors
# E2+       E4+
#"#51B8B1" "#D97D59"

## CHECK: logFC_Heatmap() looks up `cluster_levels` from the global env for
## row ordering (it's not a function argument) - loaded per-datatype below,
## right before each call, so it's always the current resolution's levels.

## fill_stat-adaptable logFC/t-stat heatmap helper (logFC_Heatmap / logFC_Heatmap_contrast)
source(here("code", "13_compile_DGE", "logFC_heatmap.R"))


#### Load DE data ####

## resolution cascade: spatial domain (Visium) -> broad cell type (sn_broad)
## -> fine subcluster (sn_fine). Xenium layers are handled in the validation section.
de_data_types <- c("Visium", "sn_broad", "sn_fine")

load_DE_data <- function(datatype) {
    DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s.Rds", datatype))
    stopifnot("DE data file not found" = file.exists(DE_data_fn))
    message("loading ", datatype, ": ", DE_data_fn)
    readRDS(DE_data_fn)
}

## each Rds already carries its own data_type column, so a plain bind is fine -
## bind_rows()/list_rbind() fill NA for columns not shared across data types
DE_data <- de_data_types |>
    map(load_DE_data) |>
    list_rbind()

message(nrow(DE_data), " rows across ", n_distinct(DE_data$data_type), " data types")
count(DE_data, data_type)

# DE_data |> filter(gene_name == "PLP1") |> arrange(vlmf_adj.P.Val)  # quick sanity check

DE_data |>
    filter(grepl("WNT", gene_name)) |>
    arrange(vlmf_adj.P.Val) |>
    select(data_type, cluster, gene_name, starts_with("vlmf"))

## quick lookup
DE_data |>
    filter(gene_name == "TRIL",
           # cluster == "vVasc",
           grepl("Astro", cluster)) |>
    select(data_type, cluster, gene_name, starts_with("vlmf")) |>
    arrange(vlmf_adj.P.Val)

#### Anchor gene panel ####

## Genes that recur across all three DE resolutions (Section 2.4-2.6) - same
## myelination/differentiation-loss + inflammatory story visible at every
## level, funneling down to concentrate in Oligo.3 at fine resolution.
## Order is fixed (not data-driven) so gene order stays identical across all
## three resolution panels for direct visual comparison.
anchor_genes <- c(
    "NPTXR",                                  # cross-cell-type signal (Astro/Oligo/Excit) at broad resolution
    "MAL", "KLK6", "CNDP1",                    # spatial-domain anchors. MAL recurs at fine resolution too
    "PLP1", "MAG", "MBP", "SOX10", "OPALIN", "LPAR1",  # Oligo.3 myelination/differentiation program, downregulated
    "FOS", "TLR2", "STAT1", "STAT4",           # Oligo.3 inflammatory/interferon program, upregulated
    "GPM6A",                                  # Oligo.3 upreg, glutamatergic-synapse
    "CPNE4",                                  # L6 multi-cell environment gene
    "FZD8", "CABLES1",                        # Mediation genes
    "TRIL",                                   # Astro.2 upreg gene with additional support
    "MAPT", "MAP2",
    "WNT7A"
)


#### Cluster colors + levels per resolution ####

## logFC_Heatmap() reads `cluster_levels` (row ordering) from wherever it's
## defined at call time - since Visium/sn_broad/sn_fine each have their own
## cluster set, this has to be reloaded per datatype rather than set once.
load_cluster_lookup <- function(datatype){
    if(datatype == "sn_broad"){
        load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
        list(cluster_colors = cell_type_colors$broad, cluster_levels = names(cell_type_colors$broad))
    } else if(datatype == "sn_fine"){
        load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
        list(cluster_colors = cell_type_colors$anno, cluster_levels = names(cell_type_colors$anno))
    } else if(datatype == "Visium"){
        load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
        list(cluster_colors = SpD_colors, cluster_levels = names(SpD_colors))
    } else if(datatype == "Xenium_cell_type_anno"){
        load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
        list(cluster_colors = cell_type_colors$anno, cluster_levels = names(cell_type_colors$anno))
    } else if(datatype == "Xenium_SpX"){
        load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
        list(cluster_colors = cell_type_colors$anno, cluster_levels = names(cell_type_colors$anno))
    } else if(datatype == "Xenium_Oligo.3_Astro"){
        stop("TODO: no cluster_colors/cluster_levels source defined yet for Xenium_Oligo.3_Astro")
    } else {
        stop(sprintf("no cluster color/level lookup defined for datatype '%s'", datatype))
    }
}


#### Build one t-stat heatmap per resolution ####

## no cluster-level filtering on the DE data itself - showing every SpD /
## broad cell type / fine subcluster at each resolution is what makes the
## funnel visible (signal diffuse at spatial-domain and broad-cell-type
## level, concentrated in Oligo.3 - and to a lesser extent Astro.3/Micro.4 -
## at fine resolution)
## using a for loop (not walk()) since cluster_colors/cluster_levels need to
## land in the global env for logFC_Heatmap() to pick them up each iteration
for(dt in de_data_types){

    lookup <- load_cluster_lookup(dt)
    cluster_colors <- lookup$cluster_colors  # logFC_Heatmap() auto-builds a row/col color annotation from this
    cluster_levels <- lookup$cluster_levels  # logFC_Heatmap() reads this from the global env

    dt_data <- DE_data |> filter(data_type == dt)

    ## full version - every cluster in this resolution
    logFC_Heatmap(
        data = dt_data,
        gene_list = anchor_genes,
        title = "anchor_genes",
        datatype = dt,
        fill_stat = "vlmf_t",
        flip = TRUE,
        order_genes = FALSE,  # keep the fixed thematic gene order above, not a data-driven one
        save = TRUE,
        h = 5.5,
        w = length(unique(dt_data$cluster))*0.5 + 2
    )

    ## filtered version - only clusters with >=1 anchor gene at adj.P.Val < 0.05
    sig_clusters <- dt_data |>
        filter(gene_name %in% anchor_genes, vlmf_adj.P.Val < 0.05) |>
        distinct(cluster) |>
        pull(cluster)

    if(length(sig_clusters) == 0){
        message(dt, ": no clusters with a significant anchor gene hit (adj.P.Val < 0.05) - skipping filtered heatmap")
    } else {
        logFC_Heatmap(
            data = dt_data |> filter(cluster %in% sig_clusters),
            gene_list = anchor_genes,
            title = "anchor_genes_sig_clusters",
            datatype = dt,
            fill_stat = "vlmf_t",
            flip = TRUE,
            order_genes = FALSE,
            save = TRUE,
            h = 5.5,
            w = length(sig_clusters) * 0.9,
            legend_side = "top"
        )
    }
}


## each call writes DGE_<datatype>_t_heatmap_anchor_genes.pdf to plot_dir
## next: pull these 3 into one composite panel (patchwork::wrap_elements +
## grid.grabExpr around each Heatmap draw call) once Xenium validation is
## added as a 4th resolution


#### Combined heatmap across all resolutions ####

## one condensed heatmap spanning all of data_types at once instead of one
## PDF per resolution - columns are clusters, restricted per-datatype to
## whichever pass `sig_filter` (default: >=1 anchor gene at adj.P.Val < 0.05;
## swap in validate-column filtering for the Xenium version below), split
## into blocks by data_type via column_split. Always the t-statistic (not
## adaptable to logFC like logFC_Heatmap() - this is a fixed combined summary
## view). No datatype prefix on column names (dropped per request) - Visium
## SpD names are simplified to just the composition tag (drop everything
## after "~"; each tag maps to exactly one domain so this doesn't collapse
## distinct domains together). Written as a function so re-running for
## Xenium/validation, or extending data_types later, is a one-line call.
combined_t_heatmap <- function(data_types,
                               de_data = DE_data,
                               gene_list = anchor_genes,
                               sig_filter = function(df) filter(df, vlmf_adj.P.Val < 0.05),
                               signif_label = function(df) mutate(df, signif = case_when(
                                   vlmf_adj.P.Val < 0.001 ~ "***",
                                   vlmf_adj.P.Val < 0.01 ~ "**",
                                   vlmf_adj.P.Val < 0.05 ~ "*",
                                   TRUE ~ "")),
                               title = "anchor_genes_sig_clusters",
                               h = 5, w_per_col = 0.5){

    ## per-datatype: figure out which clusters clear sig_filter, and their
    ## project-standard display order, before ever building the matrix - this
    ## is what lets col_order/col_split line up exactly by construction
    pieces <- map(data_types, function(dt){
        lookup <- load_cluster_lookup(dt)
        dt_data <- de_data |> filter(data_type == dt, gene_name %in% gene_list)

        if(dt == "Visium"){
            dt_data <- dt_data |> mutate(cluster = gsub("~.*", "", cluster))
            names(lookup$cluster_colors) <- gsub("~.*", "", names(lookup$cluster_colors))
            lookup$cluster_levels <- gsub("~.*", "", lookup$cluster_levels)
        }

        sig_clusters <- dt_data |> sig_filter() |> distinct(cluster) |> pull(cluster)
        cluster_order <- lookup$cluster_levels[lookup$cluster_levels %in% sig_clusters]

        list(
            data = dt_data |> filter(cluster %in% sig_clusters),
            cluster_order = cluster_order,
            cluster_colors = lookup$cluster_colors[cluster_order]
        )
    })

    combined_data <- pieces |> map("data") |> list_rbind()

    ## no datatype prefix on cluster names means this only works if they're
    ## actually unique across data_types - check rather than assume (e.g.
    ## "Vasc" is both an sn_broad cell type and the simplified Visium tag for
    ## "Vasc~Sp9D8" - would silently merge under pivot_wider if not caught)
    dup_clusters <- combined_data |> distinct(data_type, cluster) |> count(cluster) |> filter(n > 1) |> pull(cluster)
    if(length(dup_clusters) > 0){
        stop("cluster name(s) collide across data_types and would merge in pivot_wider: ",
             paste(dup_clusters, collapse = ", "))
    }

    col_order <- pieces |> map("cluster_order") |> unlist(use.names = FALSE)
    col_split <- factor(rep(data_types, times = map_int(pieces, ~length(.x$cluster_order))), levels = data_types)
    combined_cluster_colors <- pieces |> map("cluster_colors") |> unlist()

    gene_order <- gene_list[gene_list %in% combined_data$gene_name]

    t_matrix <- combined_data |>
        dplyr::select(gene_name, cluster, vlmf_t) |>
        pivot_wider(names_from = cluster, values_from = vlmf_t) |>
        column_to_rownames("gene_name") |>
        as.matrix()
    t_matrix <- t_matrix[gene_order, col_order]

    pval_matrix <- combined_data |>
        signif_label() |>
        dplyr::select(gene_name, cluster, signif) |>
        pivot_wider(names_from = cluster, values_from = signif) |>
        column_to_rownames("gene_name") |>
        as.matrix()
    pval_matrix[is.na(pval_matrix)] <- ""
    pval_matrix <- pval_matrix[gene_order, col_order]

    max_abs <- max(abs(t_matrix), na.rm = TRUE)
    my.col <- circlize::colorRamp2(
        breaks = c(-1*max_abs, 0, max_abs),
        colors = c(APOE_carrier_colors[["E2+"]], "white", APOE_carrier_colors[["E4+"]])
    )

    col_anno <- HeatmapAnnotation(cluster = col_order,
                                  col = list(cluster = combined_cluster_colors),
                                  show_legend = FALSE)

    ht <- Heatmap(t_matrix,
                  col = my.col,
                  name = "t-statistic",
                  cluster_rows = TRUE,
                  cluster_columns = FALSE,
                  column_split = col_split,          # divides the heatmap into one block per data_type
                  bottom_annotation = col_anno,      # cluster-identity color strip
                  heatmap_legend_param = list(direction = "horizontal"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(pval_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                  })

    pdf(here(plot_dir, sprintf("DGE_combined_t_heatmap_%s.pdf", title)),
        height = h, width = ncol(t_matrix) * w_per_col + 2)
    draw(ht, heatmap_legend_side = "top")
    dev.off()

    invisible(ht)
}

combined_t_heatmap(de_data_types, h = 5.5)


#### Validation data ####

## Xenium spatial validation layers. `_wSN` files already carry a `validate`
## column (direction-concordant + meets the nominal validation threshold vs.
## discovery) - use that directly rather than re-deriving significance from a
## raw p-value/FDR column.
validation_data_types <- c("Xenium_cell_type_anno", "Xenium_cell_type_anno_SpX", "Xenium_Oligo.3_Astro")

load_DE_valid_data <- function(datatype){
    DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s_wSN.Rds", datatype))
    stopifnot("DE validation data file not found" = file.exists(DE_data_fn))
    message("loading ", datatype, ": ", DE_data_fn)
    readRDS(DE_data_fn)
}

DE_valid_data <- validation_data_types |>
    map(load_DE_valid_data) |>
    list_rbind()

## MAP2 and WNT7A have no validation rows, so they drop out of validation plots
setdiff(anchor_genes, DE_valid_data$gene_name)
# [1] "MAP2"  "WNT7A"

message(nrow(DE_valid_data), " rows across ", n_distinct(DE_valid_data$data_type), " validation data types")
count(DE_valid_data, data_type)

# quick look at which clusters actually validate, per resolution
DE_valid_data |>
    filter(gene_name %in% anchor_genes, validate) |>
    group_by(data_type, cluster) |>
    summarise(n = n(), genes = paste0(unique(gene_name), collapse = ", ")) |>
    arrange(genes) |>
    print(n = 30)


#### Xenium cell type validation heatmap ####

dt <- "Xenium_cell_type_anno"
lookup <- load_cluster_lookup(dt)
cluster_colors <- lookup$cluster_colors  # logFC_Heatmap() auto-builds a row/col color annotation from this
cluster_levels <- lookup$cluster_levels  # logFC_Heatmap() reads this from the global env

dt_data <- DE_valid_data |> filter(data_type == dt)

## match the clustered order from discovery dataset (MAP2 has no validation
## rows and is silently dropped from the plot)
anchor_genes2 <- c("NPTXR","CPNE4","MAP2","TLR2","FOS","GPM6A","CABLES1","FZD8","STAT4","TRIL","STAT1","MAPT","MBP","SOX10","MAL","LPAR1","OPALIN","PLP1","MAG","KLK6","CNDP1")

## full version - every cluster in this resolution
logFC_Heatmap(
    data = dt_data,
    gene_list = anchor_genes2,
    title = "anchor_genes",
    datatype = dt,
    fill_stat = "vlmf_xenium_t",
    flip = TRUE,
    order_genes = FALSE,  # keep the fixed thematic gene order above, not a data-driven one
    save = TRUE,
    h = 5,
    w = (length(unique(dt_data$cluster))*0.3) + 2,
    valid_anno = "validate"
)

## filtered version - only clusters validation
valid_clusters <- dt_data |>
    filter(gene_name %in% anchor_genes2, validate) |>
    distinct(cluster) |>
    pull(cluster)

if(length(valid_clusters) == 0){
    message(dt, ": no validation clusters - skipping filtered heatmap")
} else {
    logFC_Heatmap(
        data = dt_data |> filter(cluster %in% valid_clusters),
        gene_list = anchor_genes2,
        title = "anchor_genes_sig_clusters",
        datatype = dt,
        fill_stat = "vlmf_xenium_t",
        flip = TRUE,
        order_genes = FALSE,
        save = TRUE,
        h = 5,
        w = length(valid_clusters) * 0.9,
        valid_anno = "validate"
    )
}


#### Xenium SpX validation heatmap ####

## SpX clusters are cell_type_anno x xSpD, keyed in the data as cluster_xSpD
## (e.g. "Oligo.3_xL5"). There is no `SpX` column - the spatial-domain column
## is `xSpD`, and its colors come from SpX_colors (names match xSpD values).
## load_cluster_lookup() has no SpX branch, so the row order and column colors
## are built here explicitly; otherwise logFC_Heatmap() would use the
## Xenium cell-type levels/colors still in the global env.
dt <- "Xenium_cell_type_anno_SpX"

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpX_colors.Rdata"), verbose = TRUE)

dt_data <- DE_valid_data |>
    filter(data_type == dt) |>
    mutate(cluster = cluster_xSpD)

## CABLES1 Oligo.3 example from the mediation analysis
dt_data |>
    filter(gene_name == "CABLES1", cell_type_anno == "Oligo.3", validate) |>
    select(gene_name, cluster, cell_type_anno, xSpD, vlmf_xenium_P.Value, vlmf_xenium_t, validate)

## one row per SpX cluster. cell_type_anno and xSpD are factors - convert to
## character before any name-based indexing (a factor index into a named
## vector uses its integer codes, which silently assigns the wrong colors).
spx_meta <- dt_data |>
    distinct(cluster, cell_type_anno, xSpD) |>
    mutate(cell_type_anno = as.character(cell_type_anno),
           xSpD = as.character(xSpD)) |>
    arrange(match(cell_type_anno, names(cell_type_colors$anno)),
            match(xSpD, names(SpX_colors))) |>
    column_to_rownames("cluster")

## display order: cell type, then spatial domain (project color order)
spx_levels <- rownames(spx_meta)

## the Xenium cell-type colors are no longer relevant here; the explicit
## column annotation below supplies the color strip
if(exists("cluster_colors")) rm(cluster_colors)
cluster_levels <- spx_levels

spx_annotation <- function(clusters){
    df <- spx_meta[clusters, c("cell_type_anno", "xSpD"), drop = FALSE]
    HeatmapAnnotation(
        df = df,
        col = list(cell_type_anno = cell_type_colors$anno[unique(df$cell_type_anno)],
                   xSpD = SpX_colors[unique(df$xSpD)])
    )
}

## full version - every cluster in this resolution
logFC_Heatmap(
    data = dt_data,
    gene_list = anchor_genes2,
    title = "anchor_genes",
    datatype = dt,
    fill_stat = "vlmf_xenium_t",
    flip = TRUE,
    col_anno = spx_annotation(spx_levels),
    order_genes = FALSE,  # keep the fixed thematic gene order above, not a data-driven one
    save = TRUE,
    h = 5,
    w = (length(spx_levels) * 0.3) + 2,
    valid_anno = "validate"
)

## filtered version - only clusters validation
valid_clusters <- dt_data |>
    filter(gene_name %in% anchor_genes2, validate) |>
    distinct(cluster) |>
    pull(cluster)
valid_clusters <- spx_levels[spx_levels %in% valid_clusters]

if(length(valid_clusters) == 0){
    message(dt, ": no validation clusters - skipping filtered heatmap")
} else {
    logFC_Heatmap(
        data = dt_data |> filter(cluster %in% valid_clusters),
        gene_list = anchor_genes2,
        title = "anchor_genes_sig_clusters",
        datatype = dt,
        col_anno = spx_annotation(valid_clusters),
        fill_stat = "vlmf_xenium_t",
        flip = TRUE,
        order_genes = FALSE,
        save = TRUE,
        h = 5,
        w = length(valid_clusters) * 0.9,
        valid_anno = "validate"
    )
}


#### Session info ####

session_info()
