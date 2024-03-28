
library("colorblindr")
library("tidyverse")
library("here")

#### cell types ####
cell_type_colors <- c(Excit = "#247FBC",
                      Inhib = "#E83E38", 
                      Oligo = "#F57A00", 
                      OPC = "#D2B037",
                      Astro = "#3BB273",
                      Micro = "#663894",
                      EndoMural = "#FF56AF",
                      Other = "#4E586A")

save(cell_type_colors, file = here("processed-data","00_project_prep","cell_type_colors.Rdata"))

#### create test plots ####
## add fake data
cell_type_levels <- names(cell_type_colors)

n = 50
n_levels <- length(cell_type_levels)
test_data <- tibble(cell_type = rep(rep(cell_type_levels, each = n),2),
                    cat = rep(c("block", "mix"), each = n*n_levels),
                    x = c(unlist(map(1:n_levels, ~rnorm(n, mean = .x))), rnorm(n*n_levels, 5, 2.5)),
                    y = c(unlist(map(1:n_levels, ~rnorm(n, mean = .x))), rnorm(n*n_levels, 5, 2.5))) |>
    mutate(cell_type = factor(cell_type, levels = cell_type_levels))

test_pallet_plots <- function(pallet, pallet_name, n=50){
    
    ## preview pallet 
    color_tab <- stack(pallet) |> 
        select(name = ind, color = values) |>
        mutate(name = factor(name, levels = names(pallet)))
    
    preview <- color_tab |>
        ggplot(aes(x= "Color", y = name , fill = name)) +
        geom_tile() +
        geom_text(aes(label = color), color = "white") +
        scale_fill_manual(values = pallet) +
        theme_classic()
    
    print(preview + geom_text(aes(label = name), color = "black", nudge_y = 0.25))
    print(cvd_grid(preview))
    
    #### pairwise ####
    color_pairwise_grid <- expand_grid(color1 = names(pallet), color2 = names(pallet)) |>
        filter(color1 != color2) |>
        mutate(num = row_number(),
               pair = map2_chr(color1, color2, ~toString(sort(c(.x, .y))))) |>
        distinct(pair, .keep_all = TRUE) |>
        mutate(pair = fct_reorder(pair, num)) |>
        pivot_longer(!c(pair, num), names_to = "x", values_to = "color") 
    
    pairwise_plot <- color_pairwise_grid |> 
        ggplot(aes(x= x, y = pair , fill = color)) +
        geom_tile() +
        scale_fill_manual(values = pallet) +
        theme_bw()
    
    print(pairwise_plot)
    print(cvd_grid(pairwise_plot))
    
    #### test data ####
    n_levels <- length(pallet)
    test_data <- tibble(class = rep(rep(names(pallet), each = n),2),
                        cat = rep(c("block", "mix"), each = n*n_levels),
                        x = c(unlist(map(1:n_levels, ~rnorm(n, mean = .x))), rnorm(n*n_levels, 5, 2.5)),
                        y = c(unlist(map(1:n_levels, ~rnorm(n, mean = .x))), rnorm(n*n_levels, 5, 2.5))) |>
        mutate(class = factor(class, levels = names(pallet)))
    
    print(scatter_test_fig <- ggplot(test_data, aes(x = x, y = y, color = class)) +
               geom_point() +
               facet_wrap(~cat) +
               scale_color_manual(values = pallet) +
               labs(title = pallet_name) +
              theme_bw()
          )
    
    print(cvd_grid(scatter_test_fig))
    
    print(density_test_fig <-test_data |>
              filter(cat == "block") |>
              ggplot(aes(x = x, fill = class)) +
              geom_density(alpha = 0.7)+
              scale_fill_manual(values = pallet)+
              labs(title = pallet_name)+
              theme_bw()
          )
    
    print(cvd_grid(density_test_fig))
    
    print(boxplot_test_fig <- test_data |>
              filter(cat == "block") |>
              ggplot(aes(x = class, y = y, fill = class)) +
              geom_boxplot()+
              scale_fill_manual(values = pallet)+
              labs(title = pallet_name)+
              theme_bw()
          )
    
    print(cvd_grid(boxplot_test_fig))
}

pdf(here("plots", "00_project_prep", "01_project_colors","cell_type_color_test_plots.pdf"))
test_pallet_plots(cell_type_colors, "Cell Colors")
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()