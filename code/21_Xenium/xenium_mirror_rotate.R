
library(here)
library(spatialLIBD)
library(ggplot2)

#' spe <- qs2::qs_read(here("code", "23_spatialLIBD_app_Xenium", "spe_xenium_app.qs2"))
#' 
#' vis_clus(spe_sub,
#'          sampleid = "Br1039",
#'          clustervar = "SpX",
#'          datatype = "Xenium",
#'          point_size = 1) +
#'     theme_bw()
#' 
#' 
#' spe_test <- purrr::map(c(90, 180, 270), ~xenium_rotate(spe, sample_id =  "Br1039", degrees = .x))
#' 
#' purrr::map2(spe_test, 
#'             c(90, 180, 270),
#'             ~vis_clus(.x,
#'                       sampleid = "Br1039",
#'                       clustervar = "SpX",
#'                       datatype = "Xenium",
#'                       point_size = 1) +
#'                 labs(subtitle = paste("rotate:", .y)) +
#'                           theme_bw()
#'             )


xenium_rotate <- function(spe, sample_id = unique(spe$sample_id)[1], degrees = 90){
    
    # 'degrees' should be a single numeric divisible by 90
    stopifnot(
        length(degrees) == 1,
        is.numeric(degrees),
        degrees %% 90 == 0
    )
    
    #   Subset
    spe_sub <- spe[, colData(spe)$sample_id == sample_id]
    
    #   Convert degrees to radians. Note that positive angles represent
    #   counter-clockwise rotations
    radians <- degrees * pi / 180
    
    #   Determine the matrix by which left-multiplication represents
    #   rotation. Then apply rotation about the origin
    rotation_mat <- matrix(
        c(
            cos(radians), sin(radians), -1 * sin(radians), cos(radians)
        ),
        nrow = 2
    )
    new_coords <- rotation_mat %*% t(spatialCoords(spe_sub))
    
    dim_max <- apply(spatialCoords(spe_sub),2, max)
    
    #   Return the "rectangle" such that its top-left corner is at the spot
    #   where its previous top-left corner was
    if (degrees %% 360 == 90) {
        new_coords <- new_coords + c(dim_max[1], 0)
    } else if (degrees %% 360 == 180) {
        new_coords <- new_coords + rev(dim_max)
    } else if (degrees %% 360 == 270) {
        new_coords <- new_coords + c(0, dim_max[2])
    }
    
    #   Ensure we have integer values for coordinates, and the correct
    #   dimnames
    new_coords <- matrix(
        as.integer(t(new_coords)),
        ncol = 2,
        dimnames = dimnames(spatialCoords(spe_sub))
    )
    
    spatialCoords(spe_sub) <- new_coords
    
    return(spe_sub)
    
}

#' spe_test_mirror <- purrr::map(c("h", "v"), ~xenium_mirror(spe, sample_id =  "Br1039", axis = .x))
#' 
#' purrr::map2(spe_test_mirror,
#'             c("h", "v"),
#'             ~vis_clus(.x,
#'                       sampleid = "Br1039",
#'                       clustervar = "SpX",
#'                       datatype = "Xenium",
#'                       point_size = 1) +
#'                 labs(subtitle = paste("mirror:", .y)) +
#'                 theme_bw()
#' )
#'
#' spe_Br1039_test <- xenium_mirror(xenium_rotate(spe, sample_id =  "Br1039", degrees = 270), axis = "v")
#' 
#' vis_clus(spe_Br1039_test,
#'          sampleid = "Br1039",
#'          clustervar = "SpX",
#'          datatype = "Xenium",
#'          point_size = 1) +
#'     theme_bw()


xenium_mirror <- function(spe, sample_id = unique(spe$sample_id)[1], axis = c("h", "v")){
    
    #   Subset
    spe_sub <- spe[, colData(spe)$sample_id == sample_id]
    
    #   Parse 'axis' into a numeric vector, since the mirror operation will
    #   be performed by multiplication
    if (axis == "h") {
        refl_vec <- c(1, -1)
    } else if (axis == "v") {
        refl_vec <- c(-1, 1)
    } else if (is.null(axis)) {
        refl_vec <- c(1, 1)
    } else {
        stop(paste("Invalid 'axis' argument:", axis))
    }
    
    #   Perform mirror
    new_coords <- refl_vec * t(spatialCoords(spe_sub))
    
    dim_max <- apply(spatialCoords(spe_sub),2, max)
 
    #   Return "rectangle" to its original location depending on the axis
    if (axis == "v") {
        new_coords <- new_coords + c(dim_max[2], 0)
    } else if (axis == "h") {
        new_coords <- new_coords + c(0, dim_max[1])
    }
    
    #   Ensure we have integer values for coordinates, and the correct
    #   dimnames
    new_coords <- matrix(
        as.integer(t(new_coords)),
        ncol = 2,
        dimnames = dimnames(spatialCoords(spe_sub))
    )
    
    #   Return a copy of the SpatialExperiment with the new coordinates
    spatialCoords(spe_sub) <- new_coords
    return(spe_sub)
    
}
