## Louise Huuki-Myers, Aug 2025
## Create plot of LDSC results

#### Set Up ####

library("tidyverse")
library("getopt")
library("here")
library("sessioninfo")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type",
      "mode", "m", "1", "character", "specificity or enrichment"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt$datatype <- "sn_fine"
# opt$datatype <- "sn_broad"
# opt$datatype <- "Visium"

# opt$mode <- "specificity"
# opt$mode <- "DEG"

plot_dir <- here("plots", "18_LDSC", "05_ldsc_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### load data ####
trait_tb <- read_csv(here("processed-data", "18_LDSC", "04_ldsc_results","gwas_trait_tb.csv"))

ldsc_data <- read_csv(here("processed-data", "18_LDSC", "04_ldsc_results", sprintf("LDSC_results-%s_%s.csv", opt$datatype, opt$mode))) |>
    mutate(z.score = Coefficient_z.score,  
           mark = ifelse(FDR < 0.05, '+', ''))

if(opt$datatype == "sn_fine"){
    
    ldsc_data <- ldsc_data |> mutate(cell_type_broad = jaffelab::ss(cluster, "\\."))
    
    ldsc_data |> count(cell_type_broad)
}

ldsc_data_filter <- ldsc_data |>
    filter(FDR <= 0.1) |> 
    group_by(cluster) |> 
    filter(any(z.score != 0)) |> 
    ungroup() |>
    group_by(gwas) |>
    filter(any(z.score != 0)) |>
    ungroup()


## dot plots
ldsc_dot <- ldsc_data |>
    ggplot(aes(x = cluster, y = gwas)) +
    geom_point(aes(color = z.score, size = -log10(p_zcore))) +
    geom_text(aes(label = mark), color = 'black', size = 5) +
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# ggsave(ldsc_dot, filename = here(plot_dir, sprintf("ldsc_dot_%s_%s.png", opt$datatype, opt$mode)), height = 6, width = 6)

ldsc_dot_filter <- ldsc_data_filter |>
    ggplot(aes(x = cluster, y = gwas)) +
    geom_point(aes(color = z.score, size = -log10(p_zcore))) +
    geom_text(aes(label = mark), color = 'black', size = 5) +
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# ggsave(ldsc_dot_filter, filename = here(plot_dir, sprintf("ldsc_dot_%s_%s_filter.png", opt$datatype, opt$mode)), height = 6, width = 6)


## tile plots
ldsc_tile <-ggplot(ldsc_data, aes(x = cluster, y = gwas, fill = z.score)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0,
                         name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label = mark), color = 'black', size = 5, vjust = 0.5) +
    ggtitle("")

# ggsave(ldsc_tile, filename = here(plot_dir, sprintf("ldsc_tile_%s_%s.png", opt$datatype, opt$mode)), height = 6, width = 6)

ldsc_filter_tile <- ggplot(ldsc_data_filter, aes(x = cluster, y = gwas, fill = z.score)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', high = 'red', midpoint = 0,
                         name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label = mark), color = 'black', size = 5, vjust = 0.5) +
    ggtitle("")

# ggsave(ldsc_filter_tile, filename = here(plot_dir, sprintf("ldsc_tile_%s_%s_filter.png", opt$datatype, opt$mode)), height = 6, width = 6)

if(opt$datatype == "sn_fine"){
    
     ldsc_dot <- ldsc_dot + facet_grid(~cell_type_broad,  scales = "free_x", space = "free_x")
     ggsave(ldsc_dot, filename = here(plot_dir, sprintf("ldsc_dot_%s_%s.png", opt$datatype, opt$mode)), height = 6, width = 12)
     
     ldsc_dot_filter <- ldsc_dot_filter + facet_grid(~cell_type_broad,  scales = "free_x", space = "free_x")
     ggsave(ldsc_dot_filter, filename = here(plot_dir, sprintf("ldsc_dot_%s_%s_filter.png", opt$datatype, opt$mode)), height = 6, width = 12)
     
     ldsc_tile <- ldsc_tile + facet_grid(~cell_type_broad,  scales = "free_x", space = "free_x")
     ggsave(ldsc_tile, filename = here(plot_dir, sprintf("ldsc_tile_%s_%s.png", opt$datatype, opt$mode)), height = 8, width = 12)
     
     ldsc_filter_tile <- ldsc_filter_tile + facet_grid(~cell_type_broad,  scales = "free_x", space = "free_x")
     ggsave(ldsc_filter_tile, filename = here(plot_dir, sprintf("ldsc_tile_%s_%s_filter.png", opt$datatype, opt$mode)), height = 8, width = 12)
     
       
} else {
    ggsave(ldsc_dot, filename = here(plot_dir, sprintf("ldsc_dot_%s_%s.png", opt$datatype, opt$mode)), height = 6, width = 6)
    ggsave(ldsc_dot_filter, filename = here(plot_dir, sprintf("ldsc_dot_%s_%s_filter.png", opt$datatype, opt$mode)), height = 6, width = 6)
    ggsave(ldsc_tile, filename = here(plot_dir, sprintf("ldsc_tile_%s_%s.png", opt$datatype, opt$mode)), height = 6, width = 6)
    ggsave(ldsc_filter_tile, filename = here(plot_dir, sprintf("ldsc_tile_%s_%s_filter.png", opt$datatype, opt$mode)), height = 6, width = 6)
    
}

