#### Louise Huuki-Myers, June 2025

### adapted from https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/48d424efba0f4cbefd08b7a2ff594cf5acd455af/code/analysis/02_visium_qc/plot_anatomical_orientation.r
### Writen by @boyiguo1
 
#' Visualize layer anatomy with three genes at once
#'
#' @param spe 
#' @param sampleid 
#'
#' @returns
#' @export
#'
#' @examples
#' 
#' 
#' library(SpatialExperiment)
#' library(ggplot2)
#' library(here)
#' 
#' plot_dir <- here("plots", "05_spe_correct_cluster", "vis_RGB_test")
#' if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
#' 
#' spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
#' rownames(spe) <- rowData(spe)$gene_name
#' 
#' test2 <- vis_RGB(spe, spatial = FALSE, point_size = 2.5)
#' ggsave(test2, filename = here(plot_dir, "test2.png"))
#' 
#' test_Br5517 <- vis_RGB(spe, sampleid = "Br5517", spatial = TRUE, point_size = 2.5)
#' ggsave(test_Br5517, filename = here(plot_dir, "test_Br5517.png"))

vis_RGB <- function(spe,
                    sampleid = unique(spe$sample_id)[1],
                    gene_red = "MBP",
                    gene_green = "PCP4",
                    gene_blue = "SNAP25",
                    spatial = TRUE,
                    assayname = "logcounts",
                    image_id = "lowres",
                    minCount = 0,
                    alpha = NA,
                    auto_crop = TRUE,
                    point_size = 2){
    
    ## Some variables
    pxl_row_in_fullres <- pxl_col_in_fullres <- key <- COUNT <- NULL
    
    spe <- spe[, spe$sample_id == sampleid]
    
    # stopifnot(all(c("pxl_col_in_fullres", "pxl_row_in_fullres", "COUNT", "key") %in% colnames(d)))
    img <-
        SpatialExperiment::imgRaster(spe, sample_id = sampleid, image_id = image_id)
    
    ## Crop the image if needed
    if (auto_crop) {
        frame_lims <-
            spatialLIBD::frame_limits(spe, sampleid = sampleid, image_id = image_id)
        img <-
            img[frame_lims$y_min:frame_lims$y_max, frame_lims$x_min:frame_lims$x_max]
        adjust <-
            list(x = frame_lims$x_min, y = frame_lims$y_min)
    } else {
        adjust <- list(x = 0, y = 0)
    }

    
    d <- as.data.frame(cbind(colData(spe), SpatialExperiment::spatialCoords(spe)), optional = TRUE)
    
    ## get normalized gene values
    geneid <- c(red = gene_red, green = gene_green, blue = gene_blue)
    stopifnot(all(geneid %in% rowData(spe)$gene_name))
    
    RGB_norm <- purrr::map(geneid, ~scales::rescale(
        logcounts(spe)[which(rowData(spe)$gene_name == .x), ] |>
            scale(center = TRUE, scale = FALSE) |>
            sapply(max, minCount), # ignore small counts - check this is correct
        to = c(0, 1)
    ))
    
    ## get RGB values
    rgb_values <- rgb(
        red = RGB_norm$red,
        green = RGB_norm$green,
        blue = RGB_norm$blue
    )
    
    #   Determine plot and legend titles
    plot_title <- paste(sampleid, "RGB")
    
    d$RGB <- rgb_values
    
    p <-
        ggplot(
            d,
            aes(
                x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id) - adjust$x,
                y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id) - adjust$y,
                fill = RGB,
                color = RGB,
                key = key
            )
        )
    
    if (spatial) {
        grob <-
            grid::rasterGrob(img,
                             width = grid::unit(1, "npc"),
                             height = grid::unit(1, "npc")
            )
        p <-
            p + spatialLIBD::geom_spatial(
                data = tibble::tibble(grob = list(grob)),
                aes(grob = grob),
                x = 0.5,
                y = 0.5
            )
    }
    
    p <- p +
        geom_point(
            shape = 21,
            size = point_size,
            stroke = 0,
            colour = "transparent",
            alpha = alpha
        ) +
        coord_fixed(expand = FALSE) +
        scale_fill_identity()
    
    
    p <- p +
        xlim(0, ncol(img)) +
        ylim(nrow(img), 0) +
        xlab("") + ylab("") +
        labs(fill = NULL, color = NULL) +
        ggtitle(plot_title) +
        theme_set(theme_bw(base_size = 20)) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(size = 10),
            legend.box.spacing = unit(0, "pt")
        )

    return(p)
    
}


vis_RGB_legend <- function(spe,
                           sampleid = unique(spe$sample_id)[1],
                           gene_red = "MBP",
                           gene_green = "PCP4",
                           gene_blue = "SNAP25",
                           assayname = "logcounts",
                           minCount = 0){
    
    spe <- spe[, spe$sample_id == sampleid]
    
    ## get normalized gene values
    geneid <- c(red = gene_red, green = gene_green, blue = gene_blue)
    stopifnot(all(geneid %in% rowData(spe)$gene_name))
    
    
    RGB_norm <- purrr::map2_dfr(geneid, names(geneid), function(gene, color){
        
        gene_tb <- tibble(gene = gene,
                          logcount = logcounts(spe)[which(rowData(spe)$gene_name == gene), ])
        
        gene_tb$scale_count <-  scales::rescale(gene_tb$logcount |>
                                                    scale(center = TRUE, scale = FALSE) |>
                                                    sapply(max, minCount), # ignore small counts
                                                to = c(0, 1)
        ) 
        
        ## get percentile     
        gene_tb_p <- gene_tb |> mutate(lc_percentile = cume_dist(logcount),
                          lc_p_round = round(lc_percentile, 1)) |>
            group_by(lc_p_round) |>
            slice_max(logcount) |>
            unique()
        
        gene_tb_p$rgb <- if(color == "red"){
            rgb(red = gene_tb_p$scale_count, green = 0, blue = 0)
        } else if(color == "green"){
            rgb(red = 0, green = gene_tb_p$scale_count, blue = 0)
        } else if(color == "blue"){
            rgb(red = 0, green = 0, blue =  gene_tb_p$scale_count)
        }
    
        return(gene_tb_p)
    })
    
    RGB_norm |> ungroup() |> count(gene)
    
    ggplot(RGB_norm, aes(x = gene, y = lc_p_round, fill = rgb)) +
        geom_tile() +
        geom_text(aes(label = round(logcount, 1)), color = "white") +
        scale_fill_identity() +
        theme_classic()
    
}

vis_RGB_legend_simple <- function(){
    
    # Create a custom legend ----
    # create random plots to get the legend
    legend_plot <- data.frame(
        x = 1, y = 1,
        Gene = c("MBP", "PCP4", "SNAP25")
    ) |>
        ggplot() +
        geom_point(aes(x, y, color = Gene)) + # Increase the size of the points
        scale_color_manual(
            name = "Gene",
            values = c(
                "MBP" = "red",
                "PCP4" = "green",
                "SNAP25" = "blue"
            ),
            guide = guide_legend(
                override.aes = list(size = 7) # Increase the size of the legend keys
            )
        ) +
        theme_void() +
        theme(
            # legend.position = "bottom",
            # legend.direction = "horizontal",
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15)
        )

    return(legend_plot)
    
}



