
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

save(cell_type_colors, file = here("processed-data","00_project_prep","cell_type_colors.Rdata"))

## APOE colors
# APOE_genotype_colors <- c(`E2/E2`="#114B5F", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46")
APOE_genotype_colors <- c(`E2/E2`="#186E8B", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46")
APOE_carrier_colors <- c(`E2+`="#398A84", `E4+`="#D46B43")

## phenotype colors
ancestry_colors <- c(EA="#1B3174",AA="#698F3F")
sex_colors <- c(M = "#5C80BC", F ="#D58BCC")

#### sample colors ####
sample_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

apoe_colors_tb <- sample_info |> 
    select(sample_id, APOE) |> 
    unique() |>
    group_by(APOE) |>
    arrange(sample_id) |>
    mutate(APOE_color = paste0(APOE, ".", row_number())) |>
    arrange(APOE_color)

apoe_colors <- DeconvoBuddies::create_cell_colors(apoe_colors_tb$APOE_color, pallet = APOE_genotype_colors, split = "\\.")
sample_colors <- apoe_colors$fine
names(sample_colors) <- apoe_colors_tb$sample_id

# Br1556    Br2582    Br5529    Br5712    Br5832    Br6161    Br5212    Br5367    Br5415    Br5426    Br5634    Br5854    Br6423
# "#186E8B" "#3E869E" "#659EB1" "#8BB6C4" "#B2CED8" "#D8E6EB" "#61C9A8" "#74CFB2" "#88D6BD" "#9CDDC8" "#B0E3D3" "#C3EADE" "#D7F1E9"
# Br6538    Br3974    Br6476    Br5276    Br5460    Br5517    Br5599    Br6085    Br6098    Br6263    Br6321    Br1039    Br1289
# "#EBF8F4" "#ED9B40" "#EEA453" "#F0AF66" "#F2B979" "#F4C38C" "#F6CD9F" "#F7D7B2" "#F9E1C5" "#FBEBD8" "#FDF5EB" "#BA3B46" "#C35760"
# Br1691    Br1706    Br2305    Br5161    Br5941
# "#CD727A" "#D78F95" "#E1ABAF" "#EBC6CA" "#F5E3E4"

## save
save(ancestry_colors, sex_colors, APOE_genotype_colors, APOE_carrier_colors, sample_colors,
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
    
    print(preview + 
              geom_text(aes(label = name), color = "black", nudge_y = 0.25) + 
              labs(title = pallet_name))
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

## plot offical colors
pdf(here("plots", "00_project_prep", "01_project_colors", "ERC_cell_type_colors.pdf"))
test_pallet_plots(cell_type_colors, "Cell Colors: ERC")
dev.off()

## other cell type colors to test
#
# ## remix for ERC
# cell_type_remix <- c(Excit = "#374194",
#                       Inhib = "#f50018", 
#                       Oligo = "#f57a00", 
#                       OPC = "#d1cf36",
#                       Astro = "#3BB273",
#                       Micro = "#7a3894",
#                       Endo = "#e8389b",
#                       Other = "#4E586A")
# 
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
# 
# cell_type_colors_list <- list(classic = cell_type_colors,
#                               remix = cell_type_remix,
#                               alt = cell_type_colors_alt,
#                               cool = cell_type_colors_cool)

# walk2(cell_type_colors_list, names(cell_type_colors_list), function(colors, name){
#     pdf(here("plots", "00_project_prep", "01_project_colors", sprintf("cell_type_color_%s_test_plots.pdf", name)))
#     test_pallet_plots(colors, paste0("Cell Colors: ", name))
#     dev.off()
# })

#### more cell colors ####
list(Crimson = "#DC143C",
     Firebrick = "#B22222",
     Indian_Red = "#CD5C5C",
     Light_Coral = "#F08080",
     Dark_Red = "#8B0000",
     Tomato = "#FF6347")

list(Red_Violet = "#C71585",
     Magenta_Red = "#D21868",
     Crimson_Red = "#DC143C",
     Pure_Red = "#FF0000",
     Scarlet_Red = "#FF4500",
     Red_Orange = "#FF6347")

list(Royal_Blue = "#4169E1",
     Dodger_Blue = "#1E90FF",
     Steel_Blue = "#4682B4",
     Sky_Blue = "#87CEEB",
     Midnight_Blue = "#191970",
     Cornflower_Blue = "#6495ED")

list(Green_Blue = "#00B3B3",
     Teal_Blue = "#0099CC",
     Cerulean_Blue = "#007ACC",
     Azure_Blue = "#005FCC",
     Cobalt_Blue = "#0040CC",
     Violet_Blue = "#2E2EFF")



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()