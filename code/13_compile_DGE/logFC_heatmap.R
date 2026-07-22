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
                          fill_stat = "vlmf_logFC"){
    
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
    
    stat_label <- switch(fill_stat,
                          "vlmf_logFC" = "log2(FC)",
                          "vlmf_t" = "t-statistic",
                          fill_stat)
    
    stat_suffix <- gsub("^vlmf_", "", fill_stat)
    
    log_fc_heatmap <- Heatmap(logFC_matrix,
                              col = my.col,
                              name = stat_label,
                              cluster_rows = FALSE,
                              cluster_columns = cluster_col,
                              right_annotation = row_anno,
                              bottom_annotation = col_anno,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                  grid.text(pval_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                              })
    
    if(save){
        pdf(here(plot_dir, sprintf("DGE_%s_%s_heatmap_%s.pdf", datatype, stat_suffix, title)), height = h, width = w)
        print(log_fc_heatmap)
        dev.off()
    } else {
        print(log_fc_heatmap)
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
                                   fill_stat = "vlmf_logFC"){
    
    if(!fill_stat %in% names(data_contrast)){
        stop(sprintf("fill_stat '%s' not a column in data_contrast", fill_stat))
    }
    
    dge_data_filter <- data_contrast |>
        filter(gene_name %in% gene_list) |>
        mutate(cluster_contrast = paste0(cluster, gsub("carrier","", contrast))) 
    
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

    log_fc_heatmap <- Heatmap(logFC_matrix,
                              col = my.col,
                              name = stat_label,
                              cluster_rows = cluster_row,
                              cluster_columns = cluster_col,
                              right_annotation = row_anno,
                              bottom_annotation = col_anno,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                  grid.text(pval_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                              })
    
    if(save){
        pdf(here(plot_dir, sprintf("DGE_%s_%s_heatmap_%s.pdf", datatype, stat_suffix, title)), height = h, width = w)
        print(log_fc_heatmap)
        dev.off()
    } else {
        print(log_fc_heatmap)
    }
    
}
