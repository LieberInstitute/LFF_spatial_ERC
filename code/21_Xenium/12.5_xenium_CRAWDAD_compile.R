## Louise Huuki-Myers, May 2026
## Compile & plot CRAWDAD results from all samples
## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/01fc3561339376fb5333163fbbcbb6f16d12ddaf/code/09_HD_cell_level/16_crawdad_plot.R

#### Set Up ####

library("here")
library("sessioninfo")
library("crawdad")
library("tidyverse")
library("spatialLIBD")

data_dir <- here("processed-data", "21_Xenium", "12.5_xenium_CRAWDAD_compile")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "12.5_xenium_CRAWDAD_compile")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Get data ####

crawdad_data_dir <- here("processed-data", "21_Xenium", "12_xenium_CRAWDAD")

min_num_signif = 5

#### custom dot plot ####

custom_dotplot = function(result_df, z_sig, filename) {
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
    
    pdf(file.path(plot_dir, filename), width = 9)
    print(p)
    dev.off()
}


