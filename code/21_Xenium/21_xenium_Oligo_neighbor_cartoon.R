## Louise Huuki-Myers, July 2026
## create example plot for xenium Oligo-Astro neighbor cartoon

#### Set Up ####

library("qs2")
library("SpatialExperiment")
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")

plot_dir <- here("plots", "21_Xenium", "21_xenium_Oligo_neighbor_cartoon")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####

message(Sys.time(), " - Load SPE data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
# filter to singlet cells
spe <- spe[,spe$spot_class == "singlet"]

table(spe$cell_type_anno == "Oligo.3")

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 

load(here("processed-data", "21_Xenium", "20_xenium_Oligo3_Astro", "Oligo3_Astro_neighbor_df.Rdata"), verbose = TRUE)

neighbor_df_details <- neighbor_df_details[neighbor_df_details$BrNum == "Br1039",]

spe <- spe[, spe$sample_id == "Br1039"]

spe_range <- function(target_barcode = "Br1039_bbfgbfki-1", radius = 50){

    target_coords <- spatialCoords(spe[, target_barcode])
    
    coord_limits <- list(x_max = target_coords[[1]] + radius,
                         x_min = target_coords[[1]] - radius,
                         y_max = target_coords[[2]] + radius,
                         y_min = target_coords[[2]] - radius
                         )
    
    
    x_in <- spatialCoords(spe)[, "x_centroid"] > coord_limits$x_min & spatialCoords(spe)[, "x_centroid"] < coord_limits$x_max
    y_in <- spatialCoords(spe)[, "y_centroid"] > coord_limits$y_min & spatialCoords(spe)[, "y_centroid"] < coord_limits$y_max
    
    spe_inrange <- spe[, x_in & y_in]
        
}

neighbor_df_details |>
    filter(neighbor_cell_type == "Astro.1", dist_class == "near",APOE_level == "high")


O3_pick <- "Br1039_bbfgbfki-1"
# O3_pick <- "Br1039_dimdldkn-1"

neighbor_df_details |> filter(reference_barcode == O3_pick)

spe_inrange <- spe_range(O3_pick, radius = 80)

single_vis_clus <- vis_clus(
    spe = spe_inrange,
    point_size = 5,
    colors = cell_type_colors$anno,
    sampleid = "Br1039",
    clustervar = "cell_type_anno",
    datatype = "Xenium"
)

ggsave(single_vis_clus, filename = here(plot_dir, "xenium_spe_range_test.png")) 


single_vis_clus <- vis_clus(
    spe = spe,
    point_size = 1.2,
    colors = cell_type_colors$anno,
    sampleid = "Br1039",
    clustervar = "cell_type_anno",
    datatype = "Xenium",
    guide_point_size = 3
) +
    geom_rect(xmin = coord_limits$x_min,
              xmax = coord_limits$x_max,
              ymin = coord_limits$y_min,
              ymax = coord_limits$y_max,
              color = "black",  
              fill = "#00000000")

ggsave(single_vis_clus, filename = here(plot_dir, "xenium_spe_Br1039_vis_cells.png"), height = 10, width = 10) 


match("APOE", colnames(colData(spe)))

colnames(colData(spe_inrange))[[17]] <- "APOE_genotype"

sqrt(spe_inrange$nucleus_area/pi)

single_vis_apoe <- vis_gene(
    spe = spe_inrange,
    sampleid = "Br1039",
    geneid = "APOE",
    point_size = 5,
    spatial = FALSE,
    datatype = "Xenium"
)

ggsave(single_vis_apoe, filename = here(plot_dir, "xenium_spe_range_APOE.png")) 

#### Plot SpX with same dims ###

spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

single_vis_clus_spx <- vis_clus(
    spe = spe,
    point_size = 1.5,
    colors = metadata(spe)$SpX_colors,
    sampleid = "Br1039",
    clustervar = "SpX",
    datatype = "Xenium",
    guide_point_size = 3
) +
    geom_rect(xmin = coord_limits$x_min,
              xmax = coord_limits$x_max,
              ymin = coord_limits$y_min,
              ymax = coord_limits$y_max,
              color = "black",  
              fill = "#00000000")

ggsave(single_vis_clus_spx, filename = here(plot_dir, "xenium_spe_Br1039_vis_SpX.png"), height = 10, width = 10) 

# slurmjobs::job_single('21_xenium_Oligo_neighbor_cartoon', create_shell = TRUE, memory = '200G', command = "Rscript 21_xenium_Oligo_neighbor_cartoon.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()



