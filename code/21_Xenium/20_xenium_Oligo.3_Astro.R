## Louise Huuki-Myers, June 2026
## Run nearest neighbors to find Oligo.3-Astro microenvironments

#### Set Up ####

library("qs2")
library("SpatialExperiment")
library("tidyverse")
library("here")
library("sessioninfo")


data_dir <- here("processed-data", "21_Xenium", "20_xenium_Oligo3_Astro")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "20_xenium_Oligo3_Astro")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### Load data ####

message(Sys.time(), " - Load SPE data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))
# filter to singlet cells
spe <- spe[,spe$spot_class == "singlet"]

table(spe$cell_type_anno == "Oligo.3")

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 

## load functions 
source(here("code", "21_Xenium", "nearest_cell_type_neighbor.R"))

#### Get Oligo.3 nearest Astro Neighbors  neighbor ####

neighbor_df <- bpmap_nearest_cell_type_neighbor(spe,
                                                reference = "Oligo.3",
                                                reference_ct_col = "cell_type_anno",
                                                neighbor = "Astro",
                                                neighbor_ct_col = "cell_type_broad",
                                                BPPARAM = BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE))


# neighbor_df <- neighbor_df |> rename(neighbor_barcode = neightbor_barcode)


pd <- colData(spe) |>
    as.data.frame() |>
    select(cell_id, cell_type_anno, cell_type_broad, SpX, BrNum, APOE_carrier)

neighbor_df_details <- neighbor_df |>
    left_join(pd |> rename(neighbor_barcode = cell_id, neighbor_cell_type = cell_type_anno, neighbor_cell_type_broad = cell_type_anno_broad))


neighbor_df_details |>
    mutate(dist_class = ifelse(distance < 50, "near", "far")) |>
    count(BrNum, dist_class) |>
    arrange(n)

#### Neighbor Identity ####

all_astro_details <- pd |>
    filter(cell_type_broad == "Astro") |>
    count(cell_type_anno) |>
    mutate(prop = n/sum(n))

neighbor_df_details |> 
    group_by(neighbor_cell_type) |>
    summarise(n = n(),
              median_dist = median(distance))

neighbor_prop <- neighbor_df_details |>
    count(BrNum, neighbor_cell_type) |>
    group_by(BrNum) |>
    mutate(prop = n/sum(n))

neighbor_boxplot <- neighbor_prop |>
    ggplot(aes(x = neighbor_cell_type, y = prop, fill =neighbor_cell_type )) +
    geom_boxplot() +
    geom_text(data = all_astro_details |> rename(neighbor_cell_type = cell_type_anno),
              aes(x = neighbor_cell_type, y = prop, label = round(prop, 2)),
              size = 3) +
    scale_fill_manual(values = cell_type_colors$anno) +
    theme_bw()

ggsave(neighbor_boxplot, filename = here(plot_dir, "Oligo3_nearest_astro_subtype_boxplot.png"))

#### Explore Oligo.3 - Astro distance ####

summary(neighbor_df_details$distance)

distance_density <- neighbor_df_details |>
    ggplot(aes(x = distance)) +
    geom_density() +
    theme_bw()

ggsave(distance_density, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_density.png"))

distance_histogram <- neighbor_df_details |>
    ggplot(aes(x = distance)) +
    geom_histogram(binwidth = 10) +
    theme_bw()

ggsave(distance_histogram, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_histogram.png"))

distance_density_subtype <- neighbor_df_details |>
    ggplot(aes(x = distance, color = neighbor_cell_type)) +
    geom_density() +
    theme_bw()

ggsave(distance_density_subtype, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_density_subtype.png"))

distance_density_SpX <- neighbor_df_details |>
    ggplot(aes(x = distance, color = SpX)) +
    geom_density() +
    theme_bw()

ggsave(distance_density_SpX, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_density_SpX.png"))

distance_boxplot_SpX <- neighbor_df_details |>
    ggplot(aes(x = SpX, y = distance, fill = SpX)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = metadata(spe)$SpX_colors) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(legend.position = "None")

ggsave(distance_boxplot_SpX, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_boxplot_SpX.png"))

distance_boxplot_donor <- neighbor_df_details |>
    ggplot(aes(x = BrNum, y = distance)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(distance_boxplot_donor, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_boxplot_donor.png"))



