compare_stats_scatter <- function(dge_tb, stat = "t", mX, mY, FDR_cut_mX = 0.2, FDR_cut_mY = 0.2, model_name){
    
    ## define vars
    statX <- paste0(mX, "_", stat)
    statY <- paste0(mY, "_", stat)
    fdrX <- paste0(mX, "_adj.P.Val")
    fdrY <- paste0(mY, "_adj.P.Val")
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "FDR<", FDR_cut_mY))
    
    # make scatter plot
    stat_scatter <- dge_tb |>
        mutate(DE_class = case_when(!!sym(fdrX) < FDR_cut_mX & !!sym(fdrY) < FDR_cut_mY ~ "sig_both",
                                    !!sym(fdrX) < FDR_cut_mX ~ paste(mX, "FDR<", FDR_cut_mX),
                                    !!sym(fdrY) < FDR_cut_mY ~ paste(mY, "FDR<", FDR_cut_mY),
                                    TRUE ~ "None")) |>
        ggplot(aes(x = !!sym(statX), y = !!sym(statY), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        labs(title = model_name, subtitle = paste(mX, "vs.", mY)) + 
        facet_wrap(~cluster) +
        theme_bw()
    
    plot_fn = sprintf("%s_%s_stat_scatter_%s_%s-v-%s.png", opt$datatype, stat, model_name, mX, mY)
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 10, width = 10)
    
    # return(t_stat_scatter)
}


compare_contrast_stats <- function(dge_tb, 
                                   stat = "t", 
                                   m = "vlmf", 
                                   FDR_cut = 0.05, 
                                   risk = FALSE, 
                                   datatype = opt$datatype,
                                   height = 10,
                                   width = 10,
                                   contrast_1,
                                   contrast_2){
    
    ## define vars
    m_stat <- paste0(m, "_", stat)
    fdr <- paste0(m, "_adj.P.Val")
    
    ## create wide data
    
    dge_tb_wide_stat <- dge_tb |>
        select(gene_id, gene_name, cluster, contrast, !!sym(m_stat)) |>
        pivot_wider(names_from = "contrast", values_from = !!sym(m_stat), names_prefix = paste0(stat, "_"))
    
    dge_tb_wide <- dge_tb |>
        select(gene_id, gene_name, cluster, contrast, !!sym(fdr)) |>
        pivot_wider(names_from = "contrast", values_from = !!sym(fdr), names_prefix = "fdr_") |> ## get wide FDR
        full_join(dge_tb_wide_stat, by = join_by(gene_id, gene_name, cluster))
    
    ## filter risk
    if(risk){
        dge_tb_wide <- dge_tb_wide |> filter(gene_name %in% AD_risk$symbol)
    }
    
    ## calc correlation
    cor <- dge_tb_wide |>
        group_by(cluster) |>
        summarise(cor = cor.test(!!sym(paste0(stat,"_",contrast_2)), !!sym(paste0(stat,"_", contrast_1)))$estimate,
                  p_val = cor.test(!!sym(paste0(stat,"_", contrast_2)), !!sym(paste0(stat,"_", contrast_1)))$p.value) |>
        ungroup() |>
        mutate(p_val_adj = p.adjust(p_val, method = "bonf"),
               signif = ifelse(p_val_adj < 0.05, "*", ""),
               anno = sprintf("cor=%.2f%s", cor, signif))
    
    
    ## create scatter plot
    stat_scatter <- compare_stats_scatter(dge_tb = dge_tb_wide, 
                                          stat = stat, 
                                          mX = contrast_1, 
                                          mY = contrast_2, 
                                          FDR_cut_mX = FDR_cut, 
                                          FDR_cut_mY = FDR_cut, 
                                          model_name = datatype) +
        geom_label(
            data = cor, ggplot2::aes(x = -Inf, y = Inf, label = anno),
            color = "black",
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        )
    
    if(risk){
        stat_scatter <- stat_scatter + geom_text_repel(aes(label = gene_name), size = 1.5)
        plot_fn = sprintf("ancestry_%s_stat_scatter_%s_risk.png", datatype, stat)
    } else {
        stat_scatter <- stat_scatter + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5)
        plot_fn = sprintf("ancestry_%s_stat_scatter_%s.png", datatype, stat)
    }
    
    
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = height, width = width)
    
    return(cor)
}


compare_stats_scatter <- function(dge_tb, stat = "t", mX, mY, FDR_cut_mX = 0.2, FDR_cut_mY = 0.2, model_name){
    
    ## define vars
    statX <- paste0(stat, "_", mX)
    statY <- paste0(stat, "_", mY)
    fdrX <- paste0("fdr_", mX)
    fdrY <- paste0("fdr_", mY)
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "FDR<", FDR_cut_mY))
    
    
    # cor <- dge_tb |>
    #     group_by(cluster) |>
    #     summarise(cor = cor(!!sym(statX), !!sym(statY))) |>
    #     mutate(anno = sprintf("cor=%.2f", cor))
    
    # make scatter plot
    dge_tb_class <- dge_tb |>
        mutate(DE_class = case_when(!!sym(fdrX) < FDR_cut_mX & !!sym(fdrY) < FDR_cut_mY ~ "sig_both",
                                    !!sym(fdrX) < FDR_cut_mX ~ paste(mX, "FDR<", FDR_cut_mX),
                                    !!sym(fdrY) < FDR_cut_mY ~ paste(mY, "FDR<", FDR_cut_mY),
                                    TRUE ~ "None")) 
    
    # return(dge_tb_class)
    
    stat_scatter <- ggplot(data = dge_tb_class, 
                           aes(x = !!sym(statX), y = !!sym(statY), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        # geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        labs(title = model_name, subtitle = paste(mX, "vs.", mY)) + 
        facet_wrap(~cluster) +
        # geom_label(
        #     data = cor, ggplot2::aes(x = -Inf, y = Inf, label = anno),
        #     color = "black",
        #     alpha = 0.5,
        #     vjust = "inward", 
        #     hjust = "inward", 
        #     size = 2.5
        # ) +
        theme_bw() +
        theme(legend.position = "bottom")
    
    return(stat_scatter)
    
}