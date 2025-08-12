
## Adapted from https://github.com/LieberInstitute/dlpfc_asd/blob/770647ed0b6c1e1dcca142af7bb2c7700581ee4a/code/06_differential_expression/13b_odds_ratio.R#L16-L74
#' de_results <- DE_data
#' gene_list <- mathys_ct_DEGs
#' 
#' gene_set_enrichment_DE(mathys_ct_DEGs,DE_data)
#' 
DE_gene_set_enrichment <- function(gene_list,
                                   de_results,
                                   fdr_cut = 0.05,
                                   test_col = "cluster",
                                   logFC_col = "vlmf_logFC",
                                   FDR_col = "vlmf_adj.P.Val") {
    
    
    test_levels <- levels(de_results[[test_col]])
    
    map_dfr(test_levels, function(test_lev){
        
        de_results_filter <- de_results |> filter(!!sym(test_col) == test_lev)
        
        ## Keep only the genes present- clean and match gene IDs
        geneList_present <- lapply(gene_list, function(x) {
            x <- x[!is.na(x)]
            x[x %in% de_results_filter$gene_name]
        })
        
        ## warn about low power for small geneLists
        geneList_length <- sapply(geneList_present, length)
        min_genes <- 25
        if (any(geneList_length < min_genes)) {
            warning(test_lev, 
                " - Gene list with n < ",
                min_genes,
                " may have low power: ",
                paste(names(geneList_length)[geneList_length < 200], collapse = ", ")
            )
        }
        
        map_dfr(c(FALSE, TRUE), function(reverse){
            ## Select direction: reverse = TRUE for downregulated
            if (reverse) {
                sig_genes <- de_results_filter[[FDR_col]] < fdr_cut & de_results_filter[[logFC_col]] < 0 #downregulated
            } else {
                sig_genes <- de_results_filter[[FDR_col]] < fdr_cut & de_results_filter[[logFC_col]] > 0 #upregulated
            }
            
            ## Fisher’s test for each gene set
            tabList <- lapply(geneList_present, function(g) {
                table(Set = factor(de_results_filter$gene_name %in% g, c(FALSE, TRUE)),
                      Sig = factor(sig_genes, c(FALSE, TRUE)))
            })
            
            enrichList <- lapply(tabList, fisher.test, alternative = "greater")
            
            enrichTab <- data.frame(
                OR = vapply(enrichList, `[[`, numeric(1), "estimate"),
                Pval = vapply(enrichList, `[[`, numeric(1), "p.value"),
                NumSig = vapply(tabList, function(x) x[2, 2], integer(1)),
                SetSize = vapply(geneList_present, length, integer(1)),
                test = paste0(test_lev, "_", DEG_regulation = ifelse(reverse, "down", "up")),
                fdr_cut = fdr_cut,
                DEG_n = sum(sig_genes),
                ID = names(geneList_present),
                stringsAsFactors = FALSE
            )
            rownames(enrichTab) <- NULL
            return(enrichTab)
        })
        
    })
    
}
