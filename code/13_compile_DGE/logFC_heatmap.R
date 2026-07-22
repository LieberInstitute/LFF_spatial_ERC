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
                          legend_side = "right",
                          signif_stat = "vlmf_adj.P.Val",
                          signif_cut = 0.05,
                          valid_anno = NULL){
    
    if(!fill_stat %in% names(data)){
        stop(sprintf("fill_stat '%s' not a column in data", fill_stat))
    }
    
    if(is.null(valid_anno)){
        if(!signif_stat %in% names(data)){
            stop(sprintf("signif_stat '%s' not a column in data", signif_stat))
        }
    } else if(!valid_anno %in% names(data)){
        stop(sprintf("valid_anno '%s' not a column in data", valid_anno))
    }
    
    logFC_matrix <- data |>
        filter(gene_name %in% gene_list) |>
        dplyr::select(cluster, gene_name, all_of(fill_stat)) |>
        pivot_wider(names_from = gene_name, values_from = all_of(fill_stat)) |>
        column_to_rownames("cluster") |>
        as.matrix()
    
    if(!is.null(valid_anno)){
        ## validation mode overrides the significance annotation entirely -
        ## an "X" wherever this boolean column is TRUE, blank elsewhere
        pval_matrix <- data |>
            filter(gene_name %in% gene_list) |>
            mutate(signif = ifelse(.data[[valid_anno]], "X", "")) |>
            dplyr::select(cluster, gene_name, signif) |>
            pivot_wider(names_from = gene_name, values_from = signif) |>
            column_to_rownames("cluster") |>
            as.matrix()
    } else {
        pval_matrix <- data |>
            filter(gene_name %in% gene_list) |>
            mutate(signif = case_when(.data[[signif_stat]] < 0.001 ~ "***",
                                      .data[[signif_stat]] < 0.01 ~ "**",
                                      .data[[signif_stat]] < signif_cut ~ "*",
                                      TRUE ~ "")
            ) |>
            dplyr::select(cluster, gene_name, signif) |>
            pivot_wider(names_from = gene_name, values_from = signif) |>
            column_to_rownames("cluster") |>
            as.matrix()
    }
    
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
    
    ## legend explaining the cell_fun text annotations - either the */**/***
    ## significance tiers, or what "X" means when valid_anno overrides them.
    ## type="grid" + graphics lets each swatch literally draw the same
    ## character used in the heatmap cells, rather than a color/point proxy.
    if(!is.null(valid_anno)){
        sig_legend <- Legend(
            title = "Validation",
            labels = valid_anno,
            type = "grid",
            graphics = list(function(x, y, w, h) grid.text("X", x, y, gp = gpar(fontsize = 10)))
        )
    } else {
        sig_legend <- Legend(
            title = "Significance",
            labels = c(paste0("p < ", signif_cut), "p < 0.01", "p < 0.001"),
            type = "grid",
            graphics = list(
                function(x, y, w, h) grid.text("*", x, y, gp = gpar(fontsize = 10)),
                function(x, y, w, h) grid.text("**", x, y, gp = gpar(fontsize = 10)),
                function(x, y, w, h) grid.text("***", x, y, gp = gpar(fontsize = 10))
            )
        )
    }
    
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
        draw(log_fc_heatmap, heatmap_legend_side = legend_side, annotation_legend_side = legend_side, annotation_legend_list = list(sig_legend))
        dev.off()
    } else {
        draw(log_fc_heatmap, heatmap_legend_side = legend_side, annotation_legend_side = legend_side, annotation_legend_list = list(sig_legend))
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
                                   legend_side = "right",
                                   signif_stat = "vlmf_adj.P.Val",
                                   signif_cut = 0.05,
                                   valid_anno = NULL){
    
    if(!fill_stat %in% names(data_contrast)){
        stop(sprintf("fill_stat '%s' not a column in data_contrast", fill_stat))
    }
    
    if(is.null(valid_anno)){
        if(!signif_stat %in% names(data_contrast)){
            stop(sprintf("signif_stat '%s' not a column in data_contrast", signif_stat))
        }
    } else if(!valid_anno %in% names(data_contrast)){
        stop(sprintf("valid_anno '%s' not a column in data_contrast", valid_anno))
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
    
    if(!is.null(valid_anno)){
        pval_matrix <- dge_data_filter |>
            mutate(signif = ifelse(.data[[valid_anno]], "X", "")) |>
            dplyr::select(cluster_contrast, gene_name, signif) |>
            pivot_wider(names_from = gene_name, values_from = signif) |>
            column_to_rownames("cluster_contrast") |>
            as.matrix()
    } else {
        pval_matrix <- dge_data_filter |>
            mutate(signif = case_when(.data[[signif_stat]] < 0.001 ~ "***",
                                      .data[[signif_stat]] < 0.01 ~ "**",
                                      .data[[signif_stat]] < signif_cut ~ "*",
                                      TRUE ~ "")
            ) |>
            dplyr::select(cluster_contrast, gene_name, signif) |>
            pivot_wider(names_from = gene_name, values_from = signif) |>
            column_to_rownames("cluster_contrast") |>
            as.matrix()
    }
    
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

    ## legend explaining the cell_fun text annotations - same convention as
    ## logFC_Heatmap()
    if(!is.null(valid_anno)){
        sig_legend <- Legend(
            title = "Validation",
            labels = valid_anno,
            type = "grid",
            graphics = list(function(x, y, w, h) grid.text("X", x, y, gp = gpar(fontsize = 10)))
        )
    } else {
        sig_legend <- Legend(
            title = "Significance",
            labels = c(paste0("p < ", signif_cut), "p < 0.01", "p < 0.001"),
            type = "grid",
            graphics = list(
                function(x, y, w, h) grid.text("*", x, y, gp = gpar(fontsize = 10)),
                function(x, y, w, h) grid.text("**", x, y, gp = gpar(fontsize = 10)),
                function(x, y, w, h) grid.text("***", x, y, gp = gpar(fontsize = 10))
            )
        )
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
        draw(log_fc_heatmap, heatmap_legend_side = legend_side, annotation_legend_side = legend_side, annotation_legend_list = list(sig_legend))
        dev.off()
    } else {
        draw(log_fc_heatmap, heatmap_legend_side = legend_side, annotation_legend_side = legend_side, annotation_legend_list = list(sig_legend))
    }
    
}
