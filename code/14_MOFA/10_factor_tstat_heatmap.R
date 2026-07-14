
#### Set up ####

library("here")
library("tidyverse")
library("sessioninfo")
library("ComplexHeatmap")

plot_dir = here("plots", "14_MOFA", "10_factor_t-stat_heatmap")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

project_colors_path = here("processed-data", "project_colors.Rdata")
tau_colors = c(`t-` = "#684F7D", `t+` = "#AFA4B6")

#### load data ####

factor_t_test_results <- read_csv(here("processed-data", "14_MOFA", "09_plot_xenium", "factor_t_test_results.csv"))

t_stat_tile <- factor_t_test_results |>
    ggplot(aes(x = covariate, y = facet_label, fill = t_stat)) +
    geom_tile() +
    geom_text(aes(label = ifelse(p_val < 0.1, round(p_val, 2), "")), color = "white") +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#D6604D", midpoint=0)+
    theme_bw() +
    labs(title = "MOFA Factor.3 Projection to Xenium")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(t_stat_tile, filename = here(plot_dir, "MOFA_factor3_xenium_t_stat_tile.png"), height = 5, width = 4)
