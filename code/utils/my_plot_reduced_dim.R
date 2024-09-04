
#' Custom Function to plot Reduced Dims
#'
#' @param spe 
#' @param dimred 
#' @param my_var 
#' @param cat_var 
#' @param save_plot 
#' @param sufix 
#' @param color_pal 
#'
#' @return
#' @export
#'
#' @examples
my_plot_reduced_dim <- function(spe,
                                dimred = "UMAP",
                                my_var = "sample_id",
                                cat_var = TRUE,
                                save_plot = TRUE,
                                sufix = NULL,
                                color_pal = NULL){
    
    rd_x = paste0(dimred, ".1")
    rd_y = paste0(dimred, ".2")
    
    rd_plot <- ggcells(spe, mapping = aes(x = !!sym(rd_x),
                                          y = !!sym(rd_y),
                                          color = !!sym(my_var)))+
        geom_point(size = 0.2, alpha = 0.3) +
        coord_equal() +
        theme_bw() +
        labs(x = paste(dimred, "Dimension 1"),
             y = paste(dimred, "Dimension 2"))
    
    if(cat_var){ ## add larger legend to catagorical variables
        rd_plot <- rd_plot + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
    } else {
        rd_plot <- rd_plot + scale_color_viridis()
    }
    
    if(!is.null(color_pal))  rd_plot <- rd_plot + scale_color_manual(values = color_pal)
    
    if(save_plot){
        plot_name <- paste(c(dimred, my_var, sufix), collapse = "_")
        plot_name <- paste0(plot_name, ".png")
        
        message("saving: ", plot_name)
        
        ggsave(rd_plot, filename = here(plot_dir, plot_name))
    } 
    
    return(rd_plot)
}