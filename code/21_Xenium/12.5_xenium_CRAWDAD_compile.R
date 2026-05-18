## Louise Huuki-Myers, May 2026
## Compile & plot CRAWDAD results from all samples
## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/01fc3561339376fb5333163fbbcbb6f16d12ddaf/code/09_HD_cell_level/16_crawdad_plot.R

#### Set Up ####

library("here")
library("sessioninfo")
library("crawdad")
library("tidyverse")
library("spatialLIBD")
library(scales)

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

