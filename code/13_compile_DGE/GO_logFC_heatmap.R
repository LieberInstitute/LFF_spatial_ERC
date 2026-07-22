

## find GO terms associated w/ parent terms
go_lookup <- function(search_term){
    
    hit <- any(map_lgl(reducedTerms_list_long, ~search_term %in% .x$parentTerm))
    
    if(!hit){
        stop(sprintf("GO parent term '%s' not in any reducedTerms_list", search_term))
    }
    
    go_terms <- map(reducedTerms_list_long, ~.x |> dplyr::filter(parentTerm == search_term) |> pull(go) |> unique())
    go_table <- map_dfr(go_terms, ~compare_clus |> filter(ID %in% .x))
    return(unique(go_table$Description))
    # return(go_table)
}

# go_terms_myelination <- go_lookup("myelination")
# go_terms_synaptic <- go_lookup("synaptic membrane")

parent_term_lookup <- function(search){
    map(reducedTerms_list_long, ~unique(.x$parentTerm[grep(search, .x$parentTerm)]))
}

parent_term_lookup("calcium")

## get genes matching GO terms
get_go_genes <- function(go_terms){
    
    go_genes <- compare_clus |>
        filter(Description %in% go_terms) |>
        select(ONTOLOGY, ID, Description, geneID) |>
        mutate(geneID = str_split(geneID, "/")) |>
        unnest_longer(geneID) |>
        unique()
    
    return(go_genes)
}

# go_genes <- get_go_genes(go_terms_myelination)

## gene can map parent term multiple times
# go_genes |> count(geneID) |> arrange(-n)
# go_genes |> count(Description) |> arrange(-n)

get_go_DE_stats <- function(go_term, contrast = FALSE, fill_stat = "vlmf_logFC"){
    
    if(!go_term %in% compare_clus$Description){
        stop(sprintf("GO term '%s' not in compare_clus", go_term))
    }
    
    if(!fill_stat %in% names(DE_data)){
        stop(sprintf("fill_stat '%s' not a column in DE_data", fill_stat))
    }
    
    go_gene <- compare_clus |>
        filter(Description == go_term) |>
        select(ONTOLOGY, ID, Description, gene_name = geneID) |>
        mutate(gene_name = str_split(gene_name, "/")) |>
        unnest_longer(gene_name) |>
        unique() 
    
    if(contrast){
        go_term_de <- go_gene |>
            left_join(DE_data |> select(cluster, contrast, gene_name, all_of(fill_stat), vlmf_adj.P.Val),
                      by = "gene_name",
                      relationship = "many-to-many") |>
            mutate(cluster = paste0(cluster, gsub("carrier","", contrast))) |>
            group_by(cluster, gene_name) |>
            arrange(vlmf_adj.P.Val) |>
            slice(1) |>
            mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                      vlmf_adj.P.Val < 0.01 ~ "**",
                                      vlmf_adj.P.Val < 0.05 ~ "*",
                                      TRUE ~ "")
            ) 
        
    } else {
        go_term_de <- go_gene |>
            left_join(DE_data |> select(cluster, gene_name, all_of(fill_stat), vlmf_adj.P.Val),
                      by = "gene_name",
                      relationship = "many-to-many") |>
            group_by(cluster, gene_name) |>
            arrange(vlmf_adj.P.Val) |>
            slice(1) |>
            mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                      vlmf_adj.P.Val < 0.01 ~ "**",
                                      vlmf_adj.P.Val < 0.05 ~ "*",
                                      TRUE ~ "")
            ) 
    }
    
    missing_col <- cluster_levels[!cluster_levels %in% go_term_de$cluster]

    
    log_fc_matrix <- go_term_de |> 
        ungroup() |>
        select(gene_name, cluster, all_of(fill_stat)) |>
        pivot_wider(names_from = "cluster", values_from = all_of(fill_stat), 
                    names_expand = TRUE) |>
        column_to_rownames("gene_name")
    
    log_fc_matrix[missing_col] <- NA
    
    log_fc_matrix <- log_fc_matrix[, cluster_levels]
    
    signif_matrix <- go_term_de |> 
        select(gene_name, cluster, signif) |>
        pivot_wider(names_from = "cluster", values_from = "signif", 
                    names_expand = TRUE)|>
        column_to_rownames("gene_name")
    
    signif_matrix[missing_col] <- NA
    signif_matrix <- signif_matrix[, cluster_levels]
    
    signif_matrix[is.na(signif_matrix)] <- ""
    
    return(list(log_fc = log_fc_matrix, 
                signif = signif_matrix,
                go_gene = go_gene,
                fill_stat = fill_stat))
    
}

# go_stats <- get_go_DE_stats("myelin assembly")
# go_stats <- get_go_DE_stats(go_term = "synaptic membrane")

# go_stats <- get_go_DE_stats("myelin assembly")
# go_stats <- get_go_DE_stats("myelination")
# get_go_DE_stats("protein autoprocessing")
# get_go_DE_stats("dendrite terminus")


get_go_DE_stats_multi <- function(go_list, contrast = FALSE, fill_stat = "vlmf_logFC"){
    
    go_stats_m <- map(go_list, ~get_go_DE_stats(.x, contrast = contrast, fill_stat = fill_stat))
    # map(go_stats_m, ~ncol(.x$log_fc))
    # drop fill_stat before transposing/rbinding - it's a repeated scalar, not a stackable matrix
    go_stats_m <- map(go_stats_m, ~.x[c("log_fc", "signif", "go_gene")])
    go_stats_m <- map(list_transpose(go_stats_m, simplify = FALSE), ~do.call("rbind",.x))
    
    go_stats_m$go_gene <-  go_stats_m$go_gene |> 
        group_by(gene_name) |>
        summarise(Description = paste0(Description, collapse = " +\n")) |>
        arrange(Description)
    
    go_stats_m$log_fc <- go_stats_m$log_fc[go_stats_m$go_gene$gene_name,]
    go_stats_m$signif <- go_stats_m$signif[go_stats_m$go_gene$gene_name,]
    go_stats_m$fill_stat <- fill_stat
    
    return(go_stats_m)
}

# go_stats_multi <- get_go_DE_stats_multi(go_list = c("myelin assembly", "protein autoprocessing", "dendrite terminus"))
# go_stats_multi <- get_go_DE_stats_multi(go_list = c("protein autoprocessing", "myelination"))

GO_logfc_Heatmap <- function(go_stats, title = NULL, cluster_rows = TRUE, fill_stat = NULL){
    
    # fall back to whatever get_go_DE_stats/_multi stashed, then to logFC for old cached objects
    if(is.null(fill_stat)) fill_stat <- go_stats$fill_stat
    if(is.null(fill_stat)) fill_stat <- "vlmf_logFC"
    
    fill_mat <- as.matrix(go_stats$log_fc)
    signif <- as.matrix(go_stats$signif)
    
    # multi_line <- any(grepl("+", go_stats$go_gene$Description))
    
    anno_table <- go_stats$go_gene |>
        # rowwise() |>
        dplyr::mutate(Description_n = ifelse(grepl("\\+", Description),
                                             Description,
                                             gsub(" ", "\n", Description))
        ) |>
        # dplyr::mutate(Description_n = Description) |>
        column_to_rownames("gene_name")
    
    anno_table <- anno_table[rownames(fill_mat),]
    
    max_abs <- max(abs(fill_mat), na.rm = TRUE)
    
    my.col <- circlize::colorRamp2(
        breaks = c(-1*max_abs,0,max_abs),
        colors = c(APOE_carrier_colors[["E2+"]], "white", APOE_carrier_colors[["E4+"]])
    )
    
    stat_label <- switch(fill_stat,
                          "vlmf_logFC" = "log2(FC)",
                          "vlmf_t" = "t-statistic",
                          fill_stat)
    
    Heatmap(fill_mat,
            col = my.col,
            name = stat_label,
            cluster_rows = cluster_rows,
            cluster_columns = FALSE,
            row_split = anno_table$Description,
            row_title_rot = 0,
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(signif[i, j], x, y, gp = gpar(fontsize = 10))
            },
            column_title = title)
    
}

parent_term_heatmap <- function(search_term, pdf_suffix = gsub(" ", "_", search_term), height = 10, width = 10, contrast = FALSE, fill_stat = "vlmf_logFC"){
    
    top_terms <- map(search_term, ~get_go_genes(go_lookup(.x))|>
                         count(Description) |>
                         arrange(-n) |>
                         head(5) |>
                         pull(Description))
    
    # TODO collapse duplicated searches
    # duplicated(top_terms)
    
    # debug
    # get_go_DE_stats_multi(top_terms[[1]], contrast = TRUE, fill_stat = fill_stat)
    
    top_terms_stats <- map(top_terms, ~get_go_DE_stats_multi(.x, contrast = contrast, fill_stat = fill_stat))
    
    # print(GO_logfc_Heatmap(top_terms_stats[[1]]))
    
    stat_suffix <- gsub("^vlmf_", "", fill_stat)
    
    if(contrast){
        filename = sprintf("GO_%s_heatmap_%s_%s-%s.pdf", stat_suffix, opt$datatype, opt$contrast, pdf_suffix)
    } else {
        filename = sprintf("GO_%s_heatmap_%s-%s.pdf", stat_suffix, opt$datatype, pdf_suffix)
    }
    
    pdf(here(plot_dir, filename), height = height, width = width)
    map2(top_terms_stats, search_term, ~print(GO_logfc_Heatmap(.x, title = .y)))
    dev.off()
}
