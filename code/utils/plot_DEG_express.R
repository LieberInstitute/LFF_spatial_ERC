#' Plot gene expression for top DE genes for one cell type
#'
#' This function plots the top n differntially expressed genes for a specified 
#' cell type.
#' The gene expression is plotted as box plot with `plot_gene_express` and adds
#' annotations to each plot.
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param stats A `data.frame()` 
#' @param cluster A `character()` target cell type to plot markers for
#' @param n_genes An `integer(1)` of number of markers you'd like to plot
#' @param rank_col The `character(1)` name of column to rank genes by in
#' `stats`.
#' @param anno_col The `character(1)` name of column containing annotation in
#' `stats`.
#' @param gene_col The `character(1)` name of column containing gene name in
#' `stats` should be the same syntax as `rownames(sce)`.
#' @param cluster_col The `character(1)` name of `colData()` column containing
#' cell type for `sce` data. It should match `cellType.target` in `stats`.
#' @param color_pal A named `character(1)` vector that contains a color palette
#' matching the `cluster` values.
#' @inheritParams plot_gene_express
#'
#' @return A `ggplot2` object created with `plot_gene_express()`. It is
#' a `scater::plotExpression()` style violin plot for selected marker genes.
#' @export
#'
#' @examples
#' ## Download the processed study data from
#'
#' ## load DE data
#' sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS"))
#' DE_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_broad", "DGE_results_carrier_sn_broad.Rds"))
#'
#' ## Plot the top markers for Astrocytes
#' plot_DEG_express_top(
#'     sce = sce_pb,
#'     stats = DE_data,
#'     cluster_col = "cell_type_broad",
#'     clus = "Astro",
#'     gene_col = "gene_name"
#' )
#' @family expression plotting functions
#' @importFrom ggplot2 ggplot geom_violin geom_text facet_wrap stat_summary
plot_DEG_express_top <- function(sce,
                                 stats,
                                 clus,
                                 n_genes = 10,
                                 pval_col = "vlmf_adj.P.Val",
                                 fc_col = "vlmf_logFC",
                                 gene_col = "gene_name",
                                 cluster_col = "cluster_broad",
                                 category_col = "APOE_carrier",
                                 mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                 color_pal = NULL,
                                 plot_points = FALSE,
                                 ncol = 2,
                                 cleanY_P = 2) {
    
    stopifnot(cluster_col %in% colnames(colData(sce)))
    stopifnot(clus %in% sce[[cluster_col]])
    stopifnot(clus %in% stats$cluster)
    
    # RCMD fix
    rank_int <- Symbol <- anno_str <- Feature <- NULL
    
    title <- paste(clus, "Top", n_genes, "DEGs")
    # message(title)
    
    max_digits <- nchar(n_genes)
    
    stopifnot(pval_col %in% colnames(stats))
    stopifnot(gene_col %in% colnames(stats))
    
    ## filter to cluster
    
    lookup <- c(
        pval_col = pval_col,
        fc_col = fc_col,
        gene_col = gene_col
    )
    
    stats_filter <- stats |>
        dplyr::rename(dplyr::all_of(lookup)) |>
        dplyr::select(cluster, gene_col, pval_col, fc_col) |>
        dplyr::filter(
            cluster == clus
        ) |>
        arrange(pval_col) |>
        mutate(
            rank_col = row_number(),
            Feature = sprintf("%s:%s",
                              stringr::str_pad(rank_col, max_digits, "left"), 
                              gene_col),
            Var1 = Feature,
            signif = case_when(pval_col < 0.001 ~"***",
                               pval_col < 0.01 ~"**",
                               pval_col < 0.05 ~"*",
                               TRUE ~ "",
            ),
            anno_str = sprintf("FDR=%.2e%s\nlogFC=%.2f", pval_col, signif, fc_col)
        ) |>
        filter(rank_col <= n_genes)
    
    # return(stats_filter)
    
    if (!any(stats_filter$gene_col %in% rownames(sce))) {
        warning("genes from gene_col don't match rownames(sce), be sure to supply the correct column from stats")
    }
    
    #### clean Y ####
    cluster_index <- sce[[cluster_col]] == clus
    sce <- sce[,cluster_index]
    
    if(is.matrix(mod)) {
        my_mod <- mod[cluster_index,]
    } else if(class(mod) == "formula"){
        my_mod <- model.matrix(mod, colData(sce))
    }
    
    assays(sce)$cleanY <- jaffelab::cleaningY(logcounts(sce),
                                              mod = my_mod,
                                              P = cleanY_P)
    
    rownames(sce) <- rowData(sce)[[gene_col]]
    
    sce <- sce[stats_filter$gene_col, ]
    rownames(sce) <- stats_filter$Feature
    
    
    pe <- plot_gene_express(
        sce = sce,
        genes = stats_filter$Feature,
        assay_name = "cleanY",
        category = category_col,
        color_pal = color_pal,
        title = title,
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
        ggplot2::geom_label(
            data = stats_filter, ggplot2::aes(x = -Inf, y = -Inf, label = anno_str),
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        )
    
    
    return(pe)
}

#' @examples
#' ## Plot the top markers for Astrocytes
#' plot_DEG_express_single(
#'     sce = sce_pb,
#'     stats = DE_data,
#'     gene = c("NXPH1", "NYAP2"),
#'     cluster_col = "cell_type_broad",
#'     clus = "Astro",
#'     gene_col = "gene_name"
#' )
#' 
plot_DEG_express <- function(sce,
                                    stats,
                                    clus = "Astro",
                                    gene = "NPTXR",
                                    pval_col = "vlmf_adj.P.Val",
                                    fc_col = "vlmf_logFC",
                                    gene_col = "gene_name",
                                    cluster_col = "registration_variable",
                                    category_col = "APOE_carrier",
                                    mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                    cleanY_P = 4,
                                    color_pal = NULL,
                                    plot_points = FALSE,
                                    ncol = 2) {
    
    stopifnot(cluster_col %in% colnames(colData(sce)))
    stopifnot(clus %in% sce[[cluster_col]])
    stopifnot(clus %in% stats$cluster)
    
    # RCMD fix
    rank_int <- Symbol <- anno_str <- NULL
    
    # title <- paste(clus, "Top", n_genes, "DEGs")
    # message(title)
    
    # max_digits <- nchar(n_genes)
    
    stopifnot(pval_col %in% colnames(stats))
    stopifnot(gene_col %in% colnames(stats))
    
    ## filter to cluster
    
    lookup <- c(
        pval_col = pval_col,
        fc_col = fc_col,
        gene_col = gene_col
    )
    
    stats_filter <- stats |>
        dplyr::rename(dplyr::all_of(lookup)) |>
        dplyr::select(cluster, gene_col, pval_col, fc_col) |>
        dplyr::filter(
            cluster == clus,
            gene_col %in% gene
        ) |>
        arrange(pval_col) |>
        mutate(
            rank_col = row_number(),
            Var1 = gene_col,
            signif = case_when(pval_col < 0.001 ~"***",
                               pval_col < 0.01 ~"**",
                               pval_col < 0.05 ~"*",
                               TRUE ~ "",
            ),
            anno_str = sprintf("FDR=%.2e%s\nlogFC=%.2f", pval_col, signif, fc_col)
        )
    
    # return(stats_filter)
    
    # if (!any(stats_filter$gene_col %in% rownames(sce))) {
    #     warning("genes from gene_col don't match rownames(sce), be sure to supply the correct column from stats")
    # }
    
    #### clean Y ####
    cluster_index <- sce[[cluster_col]] == clus
    sce <- sce[,cluster_index]
    
    if(is.matrix(mod)) {
        my_mod <- mod[cluster_index,]
    } else if(class(mod) == "formula"){
        my_mod <- model.matrix(mod, colData(sce))
    }
    
    assays(sce)$cleanY <- jaffelab::cleaningY(logcounts(sce),
                                              mod = my_mod,
                                              P = cleanY_P)
    
    ## filter data
    rownames(sce) <- rowData(sce)[[gene_col]]
    
    # sce <- sce[stats_filter$gene_col, ]
    
    pe <- plot_gene_express(
        sce = sce,
        genes = stats_filter$gene_col,
        assay_name = "cleanY",
        category = category_col,
        color_pal = color_pal,
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
        ggplot2::geom_label(
            data = stats_filter, ggplot2::aes(x = -Inf, y = -Inf, label = anno_str),
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        )
    
    
    return(pe)
}

#' @examples
#' DE_ancestry_data <- readRDS(here("processed-data", "13_compile_DGE",  "05_compile_DGE_ancestry", "sn_broad", "DGE_results_ancestry_sn_broad.Rds"))
#' sce_pb$carrier_Anc <- factor(paste(sce_pb$APOE_carrier, sce_pb$Ancestry), levels = c("E2+ AA","E4+ AA", "E2+ EA", "E4+ EA"))
#'
#' ## Plot select genes for Astrocytes
#' plot_DEG_express_contrast(
#'     sce = sce_pb,
#'     stats = DE_ancestry_data,
#'     gene = c("JUP", "RTN4RL1"),
#'     cluster_col = "cell_type_broad",
#'     clus = "Astro",
#'     gene_col = "gene_name"
#' )
#' 
plot_DEG_express_contrast <- function(sce,
                                    stats,
                                    clus = "Astro",
                                    gene = NULL,
                                    n_genes = NULL,
                                    pval_col = "vlmf_adj.P.Val",
                                    fc_col = "vlmf_logFC",
                                    gene_col = "gene_name",
                                    cluster_col = "registration_variable",
                                    category_col = "carrier_Anc",
                                    contrast_1 = "carrier_AA",
                                    contrast_2 = "carrier_EA",
                                    mod = ~0+carrier_Anc + Sex + Age + pseudo_expr_chrM_ratio,
                                    cleanY_P = 4,
                                    color_pal = NULL,
                                    plot_points = FALSE,
                                    ncol = 2,
                                    title = "") {
    
    stopifnot(cluster_col %in% colnames(colData(sce)))
    stopifnot(clus %in% sce[[cluster_col]])
    stopifnot(clus %in% stats$cluster)
    
    # RCMD fix
    rank_int <- Symbol <- anno_str <- NULL
    
    plot_title <- paste(clus, title)
    # message(title)
    
    # max_digits <- nchar(n_genes)
    
    stopifnot(pval_col %in% colnames(stats))
    stopifnot(gene_col %in% colnames(stats))
    
    ## filter to cluster
    
    lookup <- c(
        pval_col = pval_col,
        fc_col = fc_col,
        gene_col = gene_col
    )
    
    stats <- stats |>
        dplyr::rename(dplyr::all_of(lookup)) |>
        dplyr::select(cluster, contrast, gene_col, pval_col, fc_col) 
    
        stats_filter <- stats |>
            dplyr::filter(
                cluster == clus,
                gene_col %in% gene
            ) |>
            arrange(pval_col) |>
            mutate(
                rank_col = row_number(),
                Var1 = gene_col,
                signif = case_when(pval_col < 0.001 ~"***",
                                   pval_col < 0.01 ~"**",
                                   pval_col < 0.05 ~"*",
                                   TRUE ~ "",
                ),
                anno_str = sprintf("%s\nFDR=%.2e%s\nlogFC=%.2f", contrast, pval_col, signif, fc_col)
            )
        
    
    if (!any(stats_filter$gene_col %in% rownames(sce))) {
        warning("genes from gene_col don't match rownames(sce), be sure to supply the correct column from stats")
    }
    
    # return(stats_filter)
    
    #### clean Y ####
    cluster_index <- sce[[cluster_col]] == clus
    sce <- sce[,cluster_index]
    
    if(is.matrix(mod)) {
        my_mod <- mod[cluster_index,]
    } else if(class(mod) == "formula"){
        my_mod <- model.matrix(mod, colData(sce))
    }
    
    assays(sce)$cleanY <- jaffelab::cleaningY(logcounts(sce),
                                              mod = my_mod,
                                              P = cleanY_P)
    
    ## filter data
    rownames(sce) <- rowData(sce)[[gene_col]]
    
    # sce <- sce[stats_filter$gene_col, ]
    
    pe <- plot_gene_express(
        sce = sce,
        genes = unique(stats_filter$gene_col),
        assay_name = "cleanY",
        category = category_col,
        color_pal = color_pal,
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
        ggplot2::geom_label( ## contrast_1 annotations - bottom left
            data = stats_filter |> filter(contrast == contrast_1), 
            ggplot2::aes(x = -Inf, y = -Inf, label = anno_str),
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        ) +
        ggplot2::geom_label( ## contrast_1 annotations - bottom right
            data = stats_filter |> filter(contrast == contrast_2), 
            ggplot2::aes(x = Inf, y = -Inf, label = anno_str),
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        ) +
        labs(title = plot_title)
    
    
    return(pe)
}

#' ## Plot select genes for Astrocytes
#'plot_DEG_express_contrast_top(
#'    sce = sce_pb,
#'    clus = "Oligo",
#'    n_genes = 6,
#'    stats = DE_ancestry_data,
#'    cluster_col = "cell_type_broad",
#'    gene_col = "gene_name"
#')

plot_DEG_express_contrast_top <- function(sce,
                                 stats,
                                 clus = "Astro",
                                 n_genes = 5,
                                 pval_col = "vlmf_adj.P.Val",
                                 fc_col = "vlmf_logFC",
                                 gene_col = "gene_name",
                                 cluster_col = "registration_variable",
                                 category_col = "carrier_Anc",
                                 contrast_1 = "carrier_AA",
                                 contrast_2 = "carrier_EA",
                                 mod = ~0+carrier_Anc + Sex + Age + pseudo_expr_chrM_ratio,
                                 cleanY_P = 4,
                                 color_pal = NULL,
                                 plot_points = FALSE,
                                 ncol = 2) {
    
    stopifnot(cluster_col %in% colnames(colData(sce)))
    stopifnot(clus %in% sce[[cluster_col]])
    stopifnot(clus %in% stats$cluster)
    
    # RCMD fix
    rank_int <- Symbol <- anno_str <- NULL
    
    plot_title <- paste(clus)
    # message(title)
    
    max_digits <- nchar(n_genes)
    
    stats <- stats |> filter(cluster == clus)
    
    for(c in c(contrast_1, contrast_2)){
        
        top_genes <- stats |>
            filter(contrast == c) |>
            arrange(vlmf_adj.P.Val) |>
            mutate(rank = row_number(),
                   Feature = paste0(stringr::str_pad(rank, max_digits, "left"), ": ", gene_col)) |>
            filter(rank <= n_genes)
        
        print(
            plot_DEG_express_contrast(
                sce = sce_pb,
                stats = stats,
                gene = top_genes$gene_name,
                cluster_col = cluster_col,
                clus = clus,
                gene_col = gene_col,
                title = paste("-", c),
                plot_points = plot_points
            )
        )
    }
    
    
}
