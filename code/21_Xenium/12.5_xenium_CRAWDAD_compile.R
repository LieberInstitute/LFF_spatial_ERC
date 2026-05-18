## Louise Huuki-Myers, May 2026
## Compile & plot CRAWDAD results from all samples
## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/01fc3561339376fb5333163fbbcbb6f16d12ddaf/code/09_HD_cell_level/16_crawdad_plot.R

#### Set Up ####

library("here")
library("sessioninfo")
library("crawdad")
library("tidyverse")
library("spatialLIBD")
library("scales")

data_dir <- here("processed-data", "21_Xenium", "12.5_xenium_CRAWDAD_compile")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "12.5_xenium_CRAWDAD_compile")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)
cell_type_levels <- names(cell_type_colors$anno)

#### Get data ####

crawdad_data_dir <- here("processed-data", "21_Xenium", "12_xenium_CRAWDAD")

crawdad_data <- map_dfr(list.files(crawdad_data_dir, full.names = TRUE), read_csv)  |>
    mutate(
        reference = factor(reference, levels = cell_type_levels),
        neighbor = factor(neighbor, levels = cell_type_levels)
    )

#   Get the Z-score significance threshold (same in all samples)
z_sig <- crawdad_data |> filter(BrNum == "Br5460") |> correctZBonferroni()
# [1] 4.12

# crawdad_data |> 
#     filter(reference == 'Vasc.VLMC', neighbor =='Astro.4') |>
#     group_by(BrNum, neighbor, scale, reference) |>
#     summarize(Z = mean(Z)) |>
#     ggplot(aes(x = scale, y = Z, color = BrNum)) +
#     geom_point() +
#     geom_line(aes(group = BrNum))


crawdad_data |> 
    filter(reference == 'Oligo.3', neighbor =='Astro.1') |>
    group_by(id = BrNum, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |> 
    vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig)

min_num_signif = 5

crawdad_data_summary <- crawdad_data |>
    group_by(BrNum, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(BrNum, neighbor, reference) |>
    filter(scale == min(scale)) |>
    #   Retain pairs where all samples all signs of Z scores agree across
    #   samples, and significance is achieved in some sufficient number of
    #   samples 
    group_by(reference, neighbor) |>
    filter(all(Z > 0) | all(Z < 0)) |>
    filter(n() >= min_num_signif) |>
    #   Take the mean Z-score and scale across samples
    group_by(neighbor, reference) |>
    summarize(scale = mean(scale), Z = mean(Z)) |>
    ungroup() |>
    #   Cap Z-score at twice the magnitude of the significance threshold
    mutate(Z = sign(Z) * pmin(abs(Z), z_sig * 2))


crawdad_data_summary |>
    filter(reference == "Oligo.3")

crawdad_data_summary |>
    filter(reference == "Oligo.5")

crawdad_data_summary |>
    filter(neighbor == "Oligo.3")

crawdad_data_summary |>
    filter(neighbor == "Astro.4")
    

#### custom dot plot ####

custom_dotplot = function(result_df, z_sig) {
    p = ggplot(
        result_df, aes(x = reference, y = neighbor, color = Z, size = scale)
    ) +
        geom_point() +
        scale_color_gradientn(
            colors = c('blue', '#CECECE', '#CECECE', 'red'),
            values = rescale(
                c(min(result_df$Z), -1 * z_sig, z_sig, max(result_df$Z))
            )
        ) +
        scale_radius(
            trans = 'reverse',
            breaks = seq(
                min(result_df$scale), max(result_df$scale), length.out = 3
            ),
            range = c(2, 8)
        ) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    return(p)
}

dotplot <- custom_dotplot(crawdad_data_summary, z_sig)

ggsave(dotplot, filename = here(plot_dir, 'erc_xenium_CRAWDAD_dot_plot.png'), height = 12 , width = 12)


#### Visualize in Sections ####

## load spe  data 
message(Sys.time(), " - Load SPE data")
spe <- qs2::qs_read(here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2"))

spe_test <- spe[,spe$BrNum == "Br1039"]

summary(spatialCoords(spe_test)[,"x_centroid"])
summary(spatialCoords(spe_test)[,"y_centroid"])


spe_test[, spatialCoords(spe_test)[,"x_centroid"] < 850]

cell_type_vis_clus <- vis_clus(spe,
         sampleid = "Br1039",
        clustervar = "cell_type_anno",
        datatype = "Xenium",
        colors = cell_type_colors$anno, 
        point_size = 1) +
    theme_bw() +
    geom_hline(yintercept = min(spatialCoords(spe_test)[,"y_centroid"]) + c(100, 200, 500, 1000, 5000))

ggsave(cell_type_vis_clus, filename = here(plot_dir, "vis_clus_cell_type.png"), height = 10, width = 10)

plot_ref_neighbor_pair <- function(spe, sampleid = "Br1039", ref, neighbor){
    
    spe$ref_neighbor <- spe$cell_type_anno
    
    spe$ref_neighbor[!spe$ref_neighbor %in% c(ref, neighbor)] <- NA
    spe$ref_neighbor <- droplevels(spe$ref_neighbor)
    
    ref_neighbor_colors <- c("red", "blue")
    names(ref_neighbor_colors) <- c(ref, neighbor)
    
    cell_type_vis_clus <- vis_clus(spe,
                                   sampleid = sampleid,
                                   clustervar = "ref_neighbor",
                                   datatype = "Xenium",
                                   colors = ref_neighbor_colors, 
                                   point_size = 1.5)
    
    ggsave(cell_type_vis_clus, filename = here(plot_dir, sprintf("vis_clus_%s_%s_%s.png", sampleid, ref, neighbor)), height = 10, width = 10)
    
}

crawdad_data_summary |>
    filter(reference == "Oligo.3")

plot_ref_neighbor_pair(spe, ref = "Astro.1", neighbor = "Vasc.Endo")
plot_ref_neighbor_pair(spe, ref = "Oligo.3", neighbor = "Oligo.1")
plot_ref_neighbor_pair(spe, ref = "Oligo.1", neighbor = "Excit.L5.2")
plot_ref_neighbor_pair(spe, ref = "Oligo.3", neighbor = "Excit.L5.2")
