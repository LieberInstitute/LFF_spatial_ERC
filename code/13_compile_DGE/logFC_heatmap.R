logFC_Heatmap_contrast <- function(data, gene_list, title, h = 4, w = 10, cluster_col = FALSE){
    
    data <- data |>
        filter(gene_name %in% gene_list) |>
        mutate(cluster_contrast = paste0(cluster, gsub("carrier","", contrast))) 
    
    logFC_matrix <- data|>
        select(cluster_contrast, gene_name, vlmf_logFC) |>
        pivot_wider(names_from = gene_name, values_from = vlmf_logFC) |>
        column_to_rownames("cluster_contrast") |>
        as.matrix()
    
    pval_matrix <- data |>
        mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                  vlmf_adj.P.Val < 0.01 ~ "**",
                                  vlmf_adj.P.Val < 0.05 ~ "*",
                                  TRUE ~ "")
        ) |>
        select(cluster_contrast, gene_name, signif) |>
        pivot_wider(names_from = gene_name, values_from = signif) |>
        column_to_rownames("cluster_contrast") |>
        as.matrix()
    
    pval_matrix[is.na(pval_matrix)] <- ""
    
    ## reorder clusters & rows
    gene_order <- order(colMeans(logFC_matrix, na.rm = TRUE))
    
    if(all(rownames(logFC_matrix) %in% cluster_levels)){
        cluster_order <- cluster_levels[cluster_levels %in% rownames(logFC_matrix)]
    } else {
        cluster_order <- order(rownames(logFC_matrix))
    }
    
    logFC_matrix <- logFC_matrix[cluster_order, gene_order]
    pval_matrix <- pval_matrix[cluster_order, gene_order]
    
    # Heatmap(logFC_matrix, cluster_rows = FALSE, cluster_columns = FALSE)
    
    pdf(here(plot_dir, sprintf("DGE_%s_logFC_heatmap_%s.pdf", datatype, title)), height = h, width = w)
    print(Heatmap(logFC_matrix,
                  name = "log(FC)",
                  cluster_rows = FALSE,
                  cluster_columns = cluster_col,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(pval_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                  }))
    dev.off()
    
}
