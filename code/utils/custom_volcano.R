
custom_volcano <- function(data, 
                           clus, 
                           FDR_cut = 0.05, 
                           p_col = "vlmf_P.Value", 
                           fdr_col = "vlmf_adj.P.Val", 
                           lfc_col = "vlmf_logFC", 
                           model_name = "carrier",
                           save = TRUE,
                           save_size = 4,
                           text = TRUE,
                           highlight_genes = NULL){
  
  # define colors
  signif_colors <- c("purple", "blue", "red")
  names(signif_colors) <- c("both", paste("FDR<", FDR_cut) , "abs(logFC)>1" )
  
  volcano <- data |>
    filter(cluster == clus) |>
    mutate(DE_class = case_when(!!sym(fdr_col) < FDR_cut & abs(vlmf_logFC) > 1  ~ "both",
                                !!sym(fdr_col) < FDR_cut ~ paste("FDR<", FDR_cut),
                                abs(!!sym(lfc_col)) > 1 ~ "abs(logFC)>1",
                                TRUE ~ "None")) |>
    ggplot(aes(x = vlmf_logFC, y = -log10(!!sym(p_col)), color = DE_class)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_manual(values = signif_colors) +
    labs(x = "log(FC)", y = "-log10(P value)") +
    theme_bw() +
    labs(title = clus) +
    theme(legend.position = "right")
  
  if(!text){
    ## dont add text
  } else if(nrow(data) > 1000){
    volcano <- volcano + 
      # geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) 
      geom_text_repel(aes(label = ifelse(!!sym(fdr_col) < FDR_cut, gene_name, "")), size = 1.5) 
  } else {
    volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
  }
  
  if(!is.null(highlight_genes)){
    volcano <- volcano + 
      geom_text_repel(aes(label = ifelse(gene_name %in% highlight_genes, gene_name, "")), size = 1.5)
  }
  
  ## Save or return
  if(save) ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_%s_%s-%s.png", datatype, model_name, clus)), 
                  height = save_size + 0.5, 
                  width = save_size)
  else return(volcano)    
}