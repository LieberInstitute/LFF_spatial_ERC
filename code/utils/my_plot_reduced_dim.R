
#' Custom Function to plot Reduced Dims
#'
#' @param spe 
#' @param dimred 
#' @param my_var 
#' @param cat_var 
#' @param save_plot 
#' @param suffix 
#' @param color_pal 
#'
#' @return
#' @export
#'
#' @examples
my_plot_reduced_dim <- function(spe,
                                prefix, 
                                dimred = "UMAP",
                                my_var = "sample_id",
                                var_type = c("cat", "con", "express"),
                                save_plot = TRUE,
                                suffix = NULL,
                                color_pal = NULL,
                                facet = FALSE,
                                plot_dir_rd = plot_dir,
                                verbose = TRUE,
                                add_label = FALSE){
    
    rd_x = paste0(dimred, ".1")
    rd_y = paste0(dimred, ".2")
    
    rd_plot <- scater::ggcells(spe, mapping = aes(x = !!sym(rd_x),
                                          y = !!sym(rd_y),
                                          color = !!sym(my_var)))+
        geom_point(size = 0.2, alpha = 0.3) +
        coord_equal() +
        theme_bw() +
        labs(x = paste(dimred, "Dimension 1"),
             y = paste(dimred, "Dimension 2"))
    
    if(var_type == "cat"){ ## add larger legend to categorical variables
        rd_plot <- rd_plot + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
        if(facet){
            rd_plot <- rd_plot + 
                facet_wrap(as.formula(paste("~", my_var))) +
                theme(legend.position = "None") # no legend with facet
        }
    } else if(var_type == "con") {
        rd_plot <- rd_plot + viridis::scale_color_viridis()
    } else if(var_type == "express") {
        rd_plot <- rd_plot + 
            viridis::scale_color_viridis(name = "logcount") +  # label colorscheme as logcounts
            labs(title = my_var) +
            theme(plot.title = element_text(face = "italic")) # make gene name title italic
    } 
    
    if(!is.null(color_pal))  rd_plot <- rd_plot + scale_color_manual(values = color_pal)
    
    if(add_label & var_type == "cat"){
        
        rd <- as.data.frame(reducedDim(spe, dimred))
        rd_x = colnames(rd)[[1]]
        rd_y = colnames(rd)[[2]]
        
        rd$cat <- spe[[my_var]]
        
        rd_mean <- rd |> 
            group_by(cat) |> 
            # summarise(mean_x = mean(!!sym(rd_x)),
            #           mean_y = mean(!!sym(rd_y))
            #           ) |> 
            summarise(m_x = median(!!sym(rd_x)),
                      m_y = median(!!sym(rd_y))
                      )
        
        rd_plot <- rd_plot + geom_text(data = rd_mean, aes(x = m_x, y = m_y, label = cat),
                                       color = "black")
        
    }
    
    if(save_plot){
        
        if(var_type == "express"){
            plot_name <- plot_name <- paste(c(paste(c(prefix, dimred, "expression"), collapse = "_"), my_var, suffix), collapse = "-")
        } else {
            plot_name <- paste(c(paste(c(prefix, dimred), collapse = "_"), my_var, suffix), collapse = "-")
            if(facet) plot_name <- paste0(plot_name, "_facet")
        }
        
        plot_name <- paste0(plot_name, ".png")
        
        if(verbose) message(plot_name, " - " ,appendLF= FALSE)
        
        ggsave(rd_plot, filename = here(plot_dir_rd, plot_name))
    } 
    
    return(rd_plot)
}