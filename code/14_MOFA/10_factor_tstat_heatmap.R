
#### Set up ####

library("here")
library("tidyverse")
library("sessioninfo")
library("ComplexHeatmap")

plot_dir = here("plots", "14_MOFA", "10_factor_t-stat_heatmap")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

project_colors_path = here("processed-data", "project_colors.Rdata")
tau_colors = c(`t-` = "#684F7D", `t+` = "#AFA4B6")

load(here("processed-data", "SpX_colors.Rdata"), verbose = TRUE)
SpX_colors <- c(all = "black", SpX_colors)
SpX_levels <- names(SpX_colors)

subset_colors <- c(All = "black", `Astro APOE: high` = "#F89441", `Astro APOE: low` = "#7E03A8", SpX_colors)
subset_levels <- names(subset_colors)

#### load data ####

factor_t_test_results <- read_csv(here("processed-data", "14_MOFA", "09_plot_xenium", "factor_t_test_results.csv")) |>
    mutate(SpD_subset = factor(SpD_subset, levels = SpX_levels),
           enviroment = factor(case_when(
               (SpD_subset == "all" & astro_subset == "all") ~ "All",
               astro_subset == "all" ~ SpD_subset,
               SpD_subset == "all" ~ paste("Astro APOE:", astro_subset),
               TRUE ~ "ERROR"),
               levels = subset_levels),
           subset  = case_when((SpD_subset == "all" & astro_subset == "all") ~ "All",
                               astro_subset == "all" ~ "SpX",
                               SpD_subset == "all" ~ "Oligo.3 Nbr",
                               TRUE ~ "ERROR")
    )


t_stat_tile <- factor_t_test_results |>
    ggplot(aes(x = covariate, y = enviroment, fill = t_stat)) +
    geom_tile() +
    geom_text(aes(label = ifelse(p_val < 0.1, round(p_val, 2), "")), color = "white") +
    facet_wrap(~subset, ncol = 1, scales = "free_y", space = "free_y", strip.position = "right") +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#D6604D", midpoint=0)+
    theme_bw() +
    labs(title = "MOFA Factor.3", 
         subtitle = "Projection to Xenium")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave(t_stat_tile, filename = here(plot_dir, "MOFA_factor3_xenium_t_stat_tile.png"), height = 6, width = 4)
