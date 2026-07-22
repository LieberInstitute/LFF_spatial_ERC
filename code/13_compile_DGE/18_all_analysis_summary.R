## Louise Huuki-Myers July 2026
## All analysis summary
##
## REWORK: t-stat heatmaps across the DE resolution cascade (Section 2.4-2.6
## of the preprint: spatial domain -> broad cell type -> fine subcluster),
## using logFC_Heatmap() with fill_stat = "vlmf_t" for a fixed anchor gene
## panel that recurs across all three levels. Xenium spatial validation
## (Xenium_cell_type_anno / Xenium_Oligo.3_Astro / Xenium_SpX) deferred until
## the base 3-resolution figure is settled. The old ggplot tile-plot approach
## (context x gene, faceted by data_type) has been dropped in favor of one
## logFC_Heatmap() call per resolution.


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
## -> fine subcluster (sn_fine). Xenium layers deferred for now.
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

DE_data_counts <- DE_data |> 
    filter(vlmf_adj.P.Val < 0.05) |>
    group_by(gene_name) |>
    summarise(n = n(),
              n_dt = length(unique(data_type)),
              clusters = paste(unique(cluster), collapse = ",")
              ) 

DE_data_counts |>
    count(n_dt)

DE_data_counts |> filter(n_dt > 1) |> arrange(-n) |> print(n = 20)

#### Anchor gene panel ####

## Genes that recur across all three DE resolutions (Section 2.4-2.6) - same
## myelination/differentiation-loss + inflammatory story visible at every
## level, funneling down to concentrate in Oligo.3 at fine resolution.
## Order is fixed (not data-driven) so gene order stays identical across all
## three resolution panels for direct visual comparison.
anchor_genes <- c(
    "NPTXR",                                  # cross-cell-type signal (Astro/Oligo/Excit) at broad resolution
    "MAL", "KLK6",                            # spatial-domain anchors (WM.uf~Sp9D7 / Vasc~Sp9D8); MAL recurs at fine resolution too
    "PLP1", "MAG", "MBP", "SOX10", "OPALIN",  # Oligo.3 myelination/differentiation program, downregulated
    "FOS", "TLR2", "STAT1", "STAT4"           # Oligo.3 inflammatory/interferon program, upregulated
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
        # TODO: confirm cell_type_colors$fine exists alongside $broad in this Rdata - unverified assumption
        list(cluster_colors = cell_type_colors$anno, cluster_levels = names(cell_type_colors$anno))
    } else if(datatype == "Visium"){
        load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
        list(cluster_colors = SpD_colors, cluster_levels = names(SpD_colors))
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
        h = 6, w = 6          # TODO: sn_fine especially may need a taller h once you see subcluster count
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
            h = 5, w = length(sig_clusters) *0.9 ,
            legend_side = "top"
        )
    }
}


## each call writes DGE_<datatype>_t_heatmap_anchor_genes.pdf to plot_dir
## next: pull these 3 into one composite panel (patchwork::wrap_elements +
## grid.grabExpr around each Heatmap draw call) once Xenium validation is
## added as a 4th resolution


#### Combined heatmap across all resolutions ####

## one condensed heatmap spanning all of de_data_types at once instead of 3
## separate PDFs - columns are clusters, restricted per-datatype to those
## with >=1 significant anchor gene hit (same filter as the "sig_clusters"
## panels above), split into blocks by data_type via column_split. Always the
## t-statistic here (not adaptable to logFC like logFC_Heatmap() - this is a
## fixed combined summary view). Written as a function so it's a one-line
## re-run once Xenium is added as a 4th resolution: just extend de_data_types.
combined_t_heatmap <- function(data_types, gene_list = anchor_genes, title = "anchor_genes_sig_clusters", h = 5, w_per_col = 0.5){
    
    ## per-datatype: figure out which clusters clear the sig threshold, and
    ## their project-standard display order, before ever building the matrix -
    ## this is what lets col_order/col_split line up exactly by construction
    pieces <- map(data_types, function(dt){
        lookup <- load_cluster_lookup(dt)
        dt_data <- DE_data |> filter(data_type == dt)
        
        sig_clusters <- dt_data |>
            filter(gene_name %in% gene_list, vlmf_adj.P.Val < 0.05) |>
            distinct(cluster) |>
            pull(cluster)
        
        cluster_order <- lookup$cluster_levels[lookup$cluster_levels %in% sig_clusters]
        
        list(
            data = dt_data |>
                filter(cluster %in% sig_clusters, gene_name %in% gene_list) |>
                mutate(cluster_id = paste0(dt, ": ", cluster)),
            col_order = paste0(dt, ": ", cluster_order),
            n_clusters = length(cluster_order)
        )
    })
    
    combined_data <- pieces |> map("data") |> list_rbind()
    col_order <- pieces |> map("col_order") |> unlist(use.names = FALSE)
    col_split <- factor(rep(data_types, times = map_int(pieces, "n_clusters")), levels = data_types)
    
    gene_order <- gene_list[gene_list %in% combined_data$gene_name]
    
    t_matrix <- combined_data |>
        dplyr::select(gene_name, cluster_id, vlmf_t) |>
        pivot_wider(names_from = cluster_id, values_from = vlmf_t) |>
        column_to_rownames("gene_name") |>
        as.matrix()
    t_matrix <- t_matrix[gene_order, col_order]
    
    pval_matrix <- combined_data |>
        mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                  vlmf_adj.P.Val < 0.01 ~ "**",
                                  vlmf_adj.P.Val < 0.05 ~ "*",
                                  TRUE ~ "")
        ) |>
        dplyr::select(gene_name, cluster_id, signif) |>
        pivot_wider(names_from = cluster_id, values_from = signif) |>
        column_to_rownames("gene_name") |>
        as.matrix()
    pval_matrix[is.na(pval_matrix)] <- ""
    pval_matrix <- pval_matrix[gene_order, col_order]
    
    max_abs <- max(abs(t_matrix), na.rm = TRUE)
    my.col <- circlize::colorRamp2(
        breaks = c(-1*max_abs, 0, max_abs),
        colors = c(APOE_carrier_colors[["E2+"]], "white", APOE_carrier_colors[["E4+"]])
    )
    
    ht <- Heatmap(t_matrix,
                  col = my.col,
                  name = "t-statistic",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  column_split = col_split,        # divides the heatmap into one block per data_type
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

combined_t_heatmap(de_data_types)

#### Validation data ####

validation_data_types <- c("Xenium_cell_type_anno", "Xenium_SpX", "Xenium_Oligo.3_Astro")

DE_valid_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "Xenium_cell_type_anno", "DGE_results_carrier_Xenium_cell_type_anno_wSN.Rds"))

DE_valid_data |>
    filter(gene_name %in% anchor_genes,
           validate)

#### Session info ####

session_info()
