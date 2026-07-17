SpX_tstat_heatmap <- function(data,
                              gene,
                              title=gene,
                              h=6,
                              w=8,
                              cluster_col=FALSE,
                              cluster_row=FALSE,
                              datatype,
                              flip=FALSE,
                              save=TRUE,
                              order_celltypes=TRUE,
                              row_anno=NULL,
                              col_anno=NULL,
                              visium_data=NULL,
                              lookup=Sp_lookup){

    gene_data <- data|>filter(gene_name==gene)

    t_matrix <- gene_data|>
        dplyr::select(SpX,cell_type_anno,vlmf_xenium_t)|>
        pivot_wider(names_from=cell_type_anno,values_from=vlmf_xenium_t)|>
        column_to_rownames("SpX")|>
        as.matrix()

    pval_matrix <- gene_data|>
        mutate(signif=case_when(vlmf_xenium_P.Value<0.1~"X",TRUE~"")
        )|>
        dplyr::select(SpX,cell_type_anno,signif)|>
        pivot_wider(names_from=cell_type_anno,values_from=signif)|>
        column_to_rownames("SpX")|>
        as.matrix()

    pval_matrix[is.na(pval_matrix)] <- ""

    ## reorder rows (SpX) & cols (cell_type_anno)
    if(exists("SpX_levels")&&all(rownames(t_matrix)%in%SpX_levels)){
        row_order <- SpX_levels[SpX_levels%in%rownames(t_matrix)]
    } else {
        row_order <- order(rownames(t_matrix))
    }

    if(order_celltypes){
        col_order <- order(colMeans(t_matrix,na.rm=TRUE))
    } else {
        col_order <- colnames(t_matrix)
    }

    t_matrix <- t_matrix[row_order,col_order]
    pval_matrix <- pval_matrix[row_order,col_order]

    max_abs <- max(abs(t_matrix),na.rm=TRUE)

    my.col <- circlize::colorRamp2(
        breaks=c(-1*max_abs,0,max_abs),
        colors=c(APOE_carrier_colors[["E2+"]],"white",APOE_carrier_colors[["E4+"]])
    )

    if(flip){
        t_matrix <- t(t_matrix)
        pval_matrix <- t(pval_matrix)
    }

    SpX_names <- if(flip) colnames(t_matrix) else rownames(t_matrix)
    cell_types <- if(flip) rownames(t_matrix) else colnames(t_matrix)

    ## sn (discovery-side) t-stat annotation, one value per cell_type_anno, with a pch
    ## marker where vlmf_sn_adj.P.Val<0.05 - vlmf_sn_t/vlmf_sn_adj.P.Val are repeated
    ## across SpX rows for a given cell_type_anno/gene so distinct() collapses to one
    ## value per cell type
    sn_t_vec <- gene_data|>distinct(cell_type_anno,vlmf_sn_t)|>deframe()
    sn_t_vec <- sn_t_vec[cell_types]

    sn_padj_vec <- gene_data|>distinct(cell_type_anno,vlmf_sn_adj.P.Val)|>deframe()
    sn_padj_vec <- sn_padj_vec[cell_types]

    sn_max <- max(abs(sn_t_vec),na.rm=TRUE)

    sn_col <- circlize::colorRamp2(
        breaks=c(-1*sn_max,0,sn_max),
        colors=c(APOE_carrier_colors[["E2+"]],"white",APOE_carrier_colors[["E4+"]])
    )

    sn_pch <- ifelse(sn_padj_vec<0.05,8,NA)

    if(flip){
        sn_anno <- rowAnnotation(
            sn_t=anno_simple(sn_t_vec,col=sn_col,pch=sn_pch,pt_gp=gpar(col="black"))
        )
        row_anno <- if(is.null(row_anno)) sn_anno else c(row_anno,sn_anno)
    } else {
        sn_anno <- columnAnnotation(
            sn_t=anno_simple(sn_t_vec,col=sn_col,pch=sn_pch,pt_gp=gpar(col="black"))
        )
        col_anno <- if(is.null(col_anno)) sn_anno else c(col_anno,sn_anno)
    }

    ## visium (discovery-side) t-stat annotation, one value per SpX, mapped through
    ## Sp_lookup (SpX_simple -> SpD_simple -> visium cluster). SpX_simple now comes
    ## straight from the combined data frame - no more parsing SpX by regex.
    ## NOTE: assumes DE_data_visium has columns gene_name, cluster, vlmf_t,
    ## vlmf_adj.P.Val - adjust the column names below if yours differ
    if(!is.null(visium_data)){

        SpX_to_simple <- gene_data|>distinct(SpX,SpX_simple)|>deframe()
        simple_to_SpD <- lookup|>dplyr::select(SpX_simple,SpD_simple)|>deframe()
        SpD_for_rows <- simple_to_SpD[SpX_to_simple[SpX_names]]

        ## visium cluster is a compound string e.g. "Inhib~Sp09D09" - strip the
        ## "~Sp..." suffix to get the simple label that matches Sp_lookup$SpD_simple
        visium_gene <- visium_data|>
            filter(gene_name==gene)|>
            mutate(SpD_simple=gsub("~Sp.*$","",cluster))

        visium_t <- visium_gene|>dplyr::select(SpD_simple,vlmf_t)|>deframe()
        visium_t_vec <- visium_t[SpD_for_rows]
        names(visium_t_vec) <- SpX_names

        ## bail out cleanly if this gene has no usable visium match, rather than
        ## letting an all-NA vector (max -> -Inf) corrupt colorRamp2's breaks
        if(all(is.na(visium_t_vec))){

            warning(sprintf("no matching visium data for gene '%s' - skipping visium annotation",gene))

        } else {

            visium_padj <- visium_gene|>dplyr::select(SpD_simple,vlmf_adj.P.Val)|>deframe()
            visium_padj_vec <- visium_padj[SpD_for_rows]
            names(visium_padj_vec) <- SpX_names

            visium_max <- max(abs(visium_t_vec),na.rm=TRUE)

            visium_col <- circlize::colorRamp2(
                breaks=c(-1*visium_max,0,visium_max),
                colors=c(APOE_carrier_colors[["E2+"]],"white",APOE_carrier_colors[["E4+"]])
            )

            visium_pch <- ifelse(visium_padj_vec<0.05,8,NA)

            if(flip){
                visium_anno <- columnAnnotation(
                    visium_t=anno_simple(visium_t_vec,col=visium_col,pch=visium_pch,pt_gp=gpar(col="black"))
                )
                col_anno <- if(is.null(col_anno)) visium_anno else c(col_anno,visium_anno)
            } else {
                visium_anno <- rowAnnotation(
                    visium_t=anno_simple(visium_t_vec,col=visium_col,pch=visium_pch,pt_gp=gpar(col="black"))
                )
                row_anno <- if(is.null(row_anno)) visium_anno else c(row_anno,visium_anno)
            }
        }
    }

    t_heatmap <- Heatmap(t_matrix,
                        col=my.col,
                        name="t-statistic",
                        column_title=title,
                        cluster_rows=cluster_row,
                        cluster_columns=cluster_col,
                        right_annotation=row_anno,
                        bottom_annotation=col_anno,
                        cell_fun=function(j,i,x,y,width,height,fill){
                            grid.text(pval_matrix[i,j],x,y,gp=gpar(fontsize=10))
                        })

    if(save){
        pdf(here(plot_dir,sprintf("Xenium_%s_tstat_heatmap_%s.pdf",datatype,gene)),height=h,width=w)
        print(t_heatmap)
        dev.off()
    } else {
        print(t_heatmap)
    }

}
