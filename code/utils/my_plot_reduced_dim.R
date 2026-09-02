
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
                                facet = NULL,
                                plot_dir_rd = plot_dir,
                                verbose = TRUE,
                                add_label = FALSE,
                                NA_gray = FALSE,
                                void = FALSE){
    
    rd_x = paste0(dimred, ".1")
    rd_y = paste0(dimred, ".2")
    
    if(NA_gray & var_type == "express") { ## relace 0s with NA/gray in color scale 
        message("Replace 0's with NAs")
        logcounts(spe)[my_var,][logcounts(spe)[my_var,] == 0] <- NA
    } 
    
    ## create plot
    rd_plot <- scater::ggcells(spe, mapping = aes(x = !!sym(rd_x),
                                                  y = !!sym(rd_y),
                                                  color = !!sym(my_var)
    )
    )
    
    
    rd_plot <- rd_plot +
        geom_point(size = 0.2, alpha = 0.3) +
        coord_equal() +
        theme_bw() +
        labs(x = paste(dimred, "Dimension 1"),
             y = paste(dimred, "Dimension 2"))
    
    if(var_type == "cat"){ ## add larger legend to categorical variables
        rd_plot <- rd_plot + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
        
    } else if(var_type == "con") {
        rd_plot <- rd_plot + viridis::scale_color_viridis()
    } else if(var_type == "express") {
        
        # max_count <- max(logcounts(spe)[my_var,])
        # min_count <- min(logcounts(spe)[my_var,][logcounts(spe)[my_var,] > 0])
        # 
        rd_plot <- rd_plot + 
            viridis::scale_color_viridis(name = "logcount", na.value = "grey80") +  # label colorscheme as logcounts
            # scale_color_gradientn(
            #     colors = c("grey80", viridisLite::viridis(256)),
            #     values = scales::rescale(c(0, min_count, max_count)),
            #     limits = c(0, max_count)
            # ) +
            labs(title = my_var) +
            theme(plot.title = element_text(face = "italic")) # make gene name title italic
    } 
    
    if(!is.null(facet)){
        rd_plot <- rd_plot + 
            facet_wrap(as.formula(paste("~", facet))) +
            theme(legend.position = "None") # no legend with facet
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
    
    if(void){
        rd_plot <- rd_plot + theme_void()
        suffix <- paste0(suffix, "void", collapse = "_")
    }
    
    if(save_plot){
        plot_name <- plot_name <- paste(c(paste(c(prefix, dimred, var_type), collapse = "_"), my_var, suffix), collapse = "-")
        if(!is.null(facet)) plot_name <- paste0(plot_name, "_facet-", facet)
        
        plot_name <- paste0(plot_name, ".png")
        
        if(verbose) message(plot_name, " - " ,appendLF= FALSE)
        
        ggsave(rd_plot, filename = here(plot_dir_rd, plot_name))
    } 
    
    return(rd_plot)
}