logFC_Heatmap <- function(data, 
                          gene_list, 
                          title, 
                          h = 4, 
                          w = 10, 
                          cluster_col = FALSE, 
                          cluster_row = FALSE,
                          datatype,
                          flip = FALSE,
                          save = TRUE, 
                          order_genes = TRUE,
                          row_anno = NULL,
                          col_anno = NULL,
                          fill_stat = "vlmf_logFC",
                          legend_side = "right"){
    
    if(!fill_stat %in% names(data)){
        stop(sprintf("fill_stat '%s' not a column in data", fill_stat))
    }
    
    logFC_matrix <- data |>
        filter(gene_name %in% gene_list) |>
        dplyr::select(cluster, gene_name, all_of(fill_stat)) |>
        pivot_wider(names_from = gene_name, values_from = all_of(fill_stat)) |>
        column_to_rownames("cluster") |>
        as.matrix()
    
    pval_matrix <- data |>
        filter(gene_name %in% gene_list) |>
        mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                  vlmf_adj.P.Val < 0.01 ~ "**",
                                  vlmf_adj.P.Val < 0.05 ~ "*",
                                  TRUE ~ "")
        ) |>
        dplyr::select(cluster, gene_name, signif) |>
        pivot_wider(names_from = gene_name, values_from = signif) |>
        column_to_rownames("cluster") |>
        as.matrix()
    
    pval_matrix[is.na(pval_matrix)] <- ""
    
    ## reorder clusters & rows
    
    ## reorder clusters & rows
    if(order_genes){
        # message("order genes")
        gene_order <- order(colMeans(logFC_matrix, na.rm = TRUE))
    } else {
        # message("DONT order genes")
        gene_order = gene_list[gene_list %in% colnames(logFC_matrix)]
    }
    
    if(all(rownames(logFC_matrix) %in% cluster_levels)){
        cluster_order <- cluster_levels[cluster_levels %in% rownames(logFC_matrix)]
    } else {
        cluster_order <- order(rownames(logFC_matrix))
    }
    
    logFC_matrix <- logFC_matrix[cluster_order, gene_order]
    pval_matrix <- pval_matrix[cluster_order, gene_order]
    
    # Heatmap(logFC_matrix, cluster_rows = FALSE, cluster_columns = FALSE)
    
    max_abs <- max(abs(logFC_matrix), na.rm = TRUE)
    
    my.col <- circlize::colorRamp2(
        breaks = c(-1*max_abs,0,max_abs),
        colors = c(APOE_carrier_colors[["E2+"]], "white", APOE_carrier_colors[["E4+"]])
    )
    
    if(flip){
        logFC_matrix <- t(logFC_matrix)
        pval_matrix <- t(pval_matrix)
    } 
    
    # return(logFC_matrix)
    
    ## auto-build a cluster-identity color strip from `cluster_colors` (a named
    ## vector keyed by cluster name, expected in the calling env - same
    ## convention as cluster_levels) unless the caller already supplied their
    ## own row_anno/col_anno. Side depends on flip: clusters are rows normally,
    ## columns once flipped.
    if(exists("cluster_colors")){
        cluster_ids <- if(flip) colnames(logFC_matrix) else rownames(logFC_matrix)
        missing_ids <- setdiff(cluster_ids, names(cluster_colors))
        if(length(missing_ids) > 0){
            message("no cluster_colors entry for: ", paste(missing_ids, collapse = ", "), " - skipping cluster color annotation")
        } else if(flip && is.null(col_anno)){
            col_anno <- HeatmapAnnotation(cluster = cluster_ids,
                                          col = list(cluster = cluster_colors),
                                          show_legend = FALSE)
        } else if(!flip && is.null(row_anno)){
            row_anno <- rowAnnotation(cluster = cluster_ids,
                                      col = list(cluster = cluster_colors),
                                      show_legend = FALSE)
        }
    }
    
    stat_label <- switch(fill_stat,
                          "vlmf_logFC" = "log2(FC)",
                          "vlmf_t" = "t-statistic",
                          fill_stat)
    
    stat_suffix <- gsub("^vlmf_", "", fill_stat)
    
    legend_direction <- if(legend_side %in% c("top", "bottom")) "horizontal" else "vertical"
    
    log_fc_heatmap <- Heatmap(logFC_matrix,
                              col = my.col,
                              name = stat_label,
                              cluster_rows = FALSE,
                              cluster_columns = cluster_col,
                              right_annotation = row_anno,
                              bottom_annotation = col_anno,
                              heatmap_legend_param = list(direction = legend_direction),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                  grid.text(pval_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                              })
    
    if(save){
        pdf(here(plot_dir, sprintf("DGE_%s_%s_heatmap_%s.pdf", datatype, stat_suffix, title)), height = h, width = w)
        draw(log_fc_heatmap, heatmap_legend_side = legend_side)
        dev.off()
    } else {
        draw(log_fc_heatmap, heatmap_legend_side = legend_side)
    }
    
    
}



logFC_Heatmap_contrast <- function(data_contrast, 
                                   gene_list, 
                                   title, 
                                   h = 4, 
                                   w = 10, 
                                   cluster_col = FALSE,
                                   cluster_row = FALSE,
                                   flip = FALSE,
                                   save = TRUE, 
                                   order_genes = TRUE,
                                   row_anno = NULL,
                                   col_anno = NULL,
                                   fill_stat = "vlmf_logFC",
                                   legend_side = "right"){
    
    if(!fill_stat %in% names(data_contrast)){
        stop(sprintf("fill_stat '%s' not a column in data_contrast", fill_stat))
    }
    
    dge_data_filter <- data_contrast |>
        filter(gene_name %in% gene_list) |>
        mutate(cluster_contrast = paste0(cluster, gsub("carrier","", contrast))) 
    
    ## cluster_contrast (e.g. "Oligo.3E4vsE2") is what ends up as the row name -
    ## keep this lookup so cluster_colors (keyed by plain cluster name) can
    ## still be resolved per row further down
    cluster_lookup <- dge_data_filter |> distinct(cluster_contrast, cluster)
    
    logFC_matrix <- dge_data_filter|>
        dplyr::select(cluster_contrast, gene_name, all_of(fill_stat)) |>
        pivot_wider(names_from = gene_name, values_from = all_of(fill_stat)) |>
        column_to_rownames("cluster_contrast") |>
        as.matrix()
    
    pval_matrix <- dge_data_filter |>
        mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                  vlmf_adj.P.Val < 0.01 ~ "**",
                                  vlmf_adj.P.Val < 0.05 ~ "*",
                                  TRUE ~ "")
        ) |>
        dplyr::select(cluster_contrast, gene_name, signif) |>
        pivot_wider(names_from = gene_name, values_from = signif) |>
        column_to_rownames("cluster_contrast") |>
        as.matrix()
    
    pval_matrix[is.na(pval_matrix)] <- ""
    
    ## reorder clusters & rows
    if(order_genes){
        # message("order genes")
        gene_order <- order(colMeans(logFC_matrix, na.rm = TRUE))
    } else {
        # message("DONT order genes")
        gene_order = gene_list[gene_list %in% colnames(logFC_matrix)]
    }
    
    if(all(rownames(logFC_matrix) %in% cluster_levels)){
        cluster_order <- cluster_levels[cluster_levels %in% rownames(logFC_matrix)]
    } else {
        cluster_order <- order(rownames(logFC_matrix))
    }
    
    logFC_matrix <- logFC_matrix[cluster_order, gene_order]
    pval_matrix <- pval_matrix[cluster_order, gene_order]
    
    # Heatmap(logFC_matrix, cluster_rows = FALSE, cluster_columns = FALSE)
    
    max_abs <- max(abs(logFC_matrix), na.rm = TRUE)

    my.col <- circlize::colorRamp2(
        breaks = c(-1*max_abs,0,max_abs),
        colors = c(APOE_carrier_colors[["E2+"]], "white", APOE_carrier_colors[["E4+"]])
    )

    stat_label <- switch(fill_stat,
                          "vlmf_logFC" = "log2(FC)",
                          "vlmf_t" = "t-statistic",
                          fill_stat)
    
    stat_suffix <- gsub("^vlmf_", "", fill_stat)

    ## auto-build a cluster-identity color strip from `cluster_colors` (same
    ## convention as logFC_Heatmap()) unless the caller already supplied
    ## row_anno. Rows here are cluster_contrast (e.g. "Oligo.3E4vsE2"), so
    ## resolve back to plain cluster via cluster_lookup before coloring - note
    ## flip isn't applied anywhere in this function, so this stays row-side.
    if(exists("cluster_colors") && is.null(row_anno)){
        row_ids <- rownames(logFC_matrix)
        base_cluster <- cluster_lookup$cluster[match(row_ids, cluster_lookup$cluster_contrast)]
        if(anyNA(base_cluster) || any(!base_cluster %in% names(cluster_colors))){
            message("no cluster_colors entry for one or more clusters - skipping cluster color annotation")
        } else {
            row_anno <- rowAnnotation(cluster = base_cluster,
                                      col = list(cluster = cluster_colors),
                                      show_legend = FALSE)
        }
    }

    log_fc_heatmap <- Heatmap(logFC_matrix,
                              col = my.col,
                              name = stat_label,
                              cluster_rows = cluster_row,
                              cluster_columns = cluster_col,
                              right_annotation = row_anno,
                              bottom_annotation = col_anno,
                              heatmap_legend_param = list(direction = if(legend_side %in% c("top", "bottom")) "horizontal" else "vertical"),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                  grid.text(pval_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                              })
    
    if(save){
        pdf(here(plot_dir, sprintf("DGE_%s_%s_heatmap_%s.pdf", datatype, stat_suffix, title)), height = h, width = w)
        draw(log_fc_heatmap, heatmap_legend_side = legend_side)
        dev.off()
    } else {
        draw(log_fc_heatmap, heatmap_legend_side = legend_side)
    }
    
}
