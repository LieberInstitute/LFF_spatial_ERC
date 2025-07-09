

sce_pseudobulk <- function(sce, cluster){
   
   ## add syntactic APOE vars
    sce$APOE_syn <- factor(gsub("/", ".", sce$APOE))
    levels(sce$APOE_syn)
    
    sce$APOE_carrier_syn <- factor(gsub("\\+", "", sce$APOE_carrier))
    levels(sce$APOE_carrier_syn)
    
    ## add E4/E4 variable
    sce$APOE_E4E4 <- sce$APOE == "E4/E4"
    table(sce$APOE_E4E4)
    
    ## get mito genes
    mito_genes <- as.logical(seqnames(sce) == "chrM")
    table(mito_genes)
    
    table(sce[[cluster]])
    
    #### run pseudobulk - cell_type_anno ####
    message(Sys.time(), " PSEUDOBULK: ", cluster)
    sce_pseudo <- registration_pseudobulk(
        sce,
        var_registration = cluster,
        var_sample_id = "sample_id",
        covars = NULL,
        min_ncells = 10,
        pseudobulk_rds_file = NULL,
        filter_expr = FALSE,
        mito_gene = mito_genes
    )
    
    message(Sys.time(), " - Done pseudobulk")
    
    message(sprintf("nrow: %d, ncol: %d", nrow(sce_pseudo), ncol(sce_pseudo)))
    
    #### Check n samples for each cell type ####
    table(sce_pseudo$APOE_carrier, sce_pseudo$cell_type_broad)
    table(sce_pseudo$APOE_carrier, sce_pseudo$registration_variable)
    
    cell_type_count <- colData(sce_pseudo) |> as.data.frame() |> dplyr::count(APOE_carrier, registration_variable)
    
    ## cell type must have two or more samples on either side of DEG split (APOE carrier)
    enough_samples <- cell_type_count |>
        filter(n >=2) |>
        dplyr::count(registration_variable) |>
        filter(n >=2) |>
        dplyr::pull(registration_variable)
    
    message("Too few samples in: ", 
            paste(levels(sce_pseudo$registration_variable)[!levels(sce_pseudo$registration_variable) %in% enough_samples], collapse = ", ")
    )
    
    ## drop too few sample cell types
    sce_pseudo <- sce_pseudo[, sce_pseudo$registration_variable %in% enough_samples]
    sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)
    
    #### Additional edits + Save ####
    ## drop all NA cols
    all_na <- sapply(colData(sce_pseudo), function(x)all(is.na(x)))
    colData(sce_pseudo) <- colData(sce_pseudo)[, names(all_na)[!all_na]]
    
    return(sce_pseudo)
    
}