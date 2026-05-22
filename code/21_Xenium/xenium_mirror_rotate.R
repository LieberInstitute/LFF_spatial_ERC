
library(here)
library(spatialLIBD)
library(ggplot2)

spe <- qs2::qs_read(here("code", "23_spatialLIBD_app_Xenium", "spe_xenium_app.qs2"))

vis_clus(spe_sub,
         sampleid = "Br1039",
         clustervar = "SpX",
         datatype = "Xenium",
         point_size = 1) +
    theme_bw()


spe_test <- purrr::map(c(90, 180, 270), ~xenium_rotate(spe, sample_id =  "Br1039", degrees = .x))

purrr::map2(spe_test, 
            c(90, 180, 270),
            ~vis_clus(.x,
                      sampleid = "Br1039",
                      clustervar = "SpX",
                      datatype = "Xenium",
                      point_size = 1) +
                labs(subtitle = paste("rotate:", .y)) +
                          theme_bw()
            )


xenium_rotate <- function(spe, sample_id, degrees = 90){
    
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
        new_coords <- new_coords + c(dim_x_max[1], 0)
    } else if (degrees %% 360 == 180) {
        new_coords <- new_coords + rev(dim_max)
    } else if (degrees %% 360 == 270) {
        new_coords <- new_coords + c(0, dim_max[2])
    }
    
    #   Ensure we have integer values for coordinates, and the correct
    #   dimnames
    new_coords <- matrix(
        as.integer(round(t(new_coords))),
        ncol = 2,
        dimnames = dimnames(spatialCoords(spe_sub))
    )
    
    spatialCoords(spe_sub) <- new_coords
    
    return(spe_sub)
    
}