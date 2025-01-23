
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
                      Endo = "#FF56AF",
                      Other = "#4E586A")

# cell_type_colors_alt <- c(Excit = "#246EBD",
#                       Inhib = "#E83870",
#                       Oligo = "#DD900B",
#                       OPC = "#BAD136",
#                       Astro = "#3BB391",
#                       Micro = "#4C3894",
#                       Endo = "#FA63C8",
#                       Other = "#4E586A")
# 
# cell_type_colors_cool <- c(Excit = "#2423BC",
#                            Inhib = "#94378A",
#                            Oligo = "#E8A738",
#                            OPC = "#DDF500",
#                            Astro = "#8CC82D",
#                            Micro = "#3BA9B2",
#                            Endo = "#5C537F")

save(cell_type_colors, file = here("processed-data","00_project_prep","cell_type_colors.Rdata"))

## APOE colors
# APOE_genotype_colors <- c(`E2/E2`="#114B5F", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46")
APOE_genotype_colors <- c(`E2/E2`="#186E8B", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46")
APOE_carrier_colors <- c(`E2+`="#398A84", `E4+`="#D46B43")

## phenotype colors
ancestry_colors <- c(EA="#1B3174",AA="#698F3F")
sex_colors <- c(M = "#5C80BC", F ="#D58BCC")

## save
save(ancestry_colors, sex_colors, APOE_genotype_colors, APOE_carrier_colors,
     file = here("processed-data", "project_colors.Rdata"))

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