## Louise Huuki-Myers, June 2026
## Run nearest neighbors to find Oligo.3-Astro microenvironments

#### Set Up ####

library("qs2")
library("SpatialExperiment")
library("tidyverse")
library("here")
library("sessioninfo")
library("DeconvoBuddies")

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

#### Build neighbor detail table ####
pd <- colData(spe) |>
    as.data.frame() |>
    select(cell_id, cell_type_anno, cell_type_broad, SpX, BrNum, APOE_carrier)

(median_dist <- median(neighbor_df$distance)) # median = 40.8

neighbor_df$Astro_APOE <- logcounts(spe)["APOE", neighbor_df_details$neighbor_barcode]

neighbor_df_details <- neighbor_df |>
    left_join(pd |> select(neighbor_barcode = cell_id, 
                           neighbor_cell_type = cell_type_anno, 
                           neighbor_cell_type_broad = cell_type_broad,
                           APOE_carrier)) |>
    left_join(pd |> select(reference_barcode = cell_id, 
                           reference_SpX = SpX))  |>
    mutate(dist_class = ifelse(distance < median_dist, "near", "far"))

neighbor_APOE_cutoff <- neighbor_df_details |>
    select(neighbor_barcode, neighbor_cell_type, Astro_APOE) |>
    unique() |> ## rm repeats
    group_by(neighbor_cell_type) |>
    summarise(median_APOE = median(Astro_APOE))

## annotate APOE expression
neighbor_df_details <- neighbor_df_details |>
    left_join(neighbor_APOE_cutoff) |>
    mutate(APOE_level = ifelse(Astro_APOE < median_APOE, "low", "high"))

neighbor_df_details |> count(APOE_level, dist_class)

neighbor_df_details |> count(BrNum, APOE_level, dist_class, reference_SpX) |> count(APOE_level, dist_class, reference_SpX, n < 10)

neighbor_df_details |> 
    count(BrNum, APOE_carrier, APOE_level, dist_class, reference_SpX) |> 
    group_by(APOE_carrier, APOE_level, dist_class, reference_SpX) |>
    summarise(n10 = sum(n > 10), 
              min_cells = min(n), 
              median_cells = median(n), 
              max_cells = max(n))


save(neighbor_df_details, file = here(data_dir, "Oligo3_Astro_neighbor_df.Rdata"))

neighbor_df_details |>
    count(BrNum, dist_class, APOE_level) |>
    arrange(-n)

#### Summarize by SpX ####

neighbor_detail_summary <- neighbor_df_details |>
    mutate(Oligo.3_Astro = paste0("nnA_", dist_class, "_APOE_", APOE_level)) |>
    group_by(Oligo.3_Astro, dist_class, APOE_level, reference_SpX) |>
    summarize(n = n()) |>
    group_by(reference_SpX) |>
    mutate(SpX_prop = n/sum(n))|>
    group_by(Oligo.3_Astro) |>
    mutate(type_prop = n/sum(n))

write.csv(neighbor_detail_summary, file = here(data_dir, "Oligo3_Astro_neighbor_SpX_summary.csv"))


type_SpX_n_tile <- neighbor_detail_summary |>
    ggplot(aes(x = Oligo.3_Astro, y = reference_SpX, fill = n)) +
    geom_tile() +
    scale_fill_gradientn(colors = c(viridisLite::plasma(100))) +
    theme_bw() +
    scale_y_discrete(limits=rev) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(type_SpX_n_tile, filename = here(plot_dir, "Oligo3_nearest_astro_type_SpX_n_tile.png"))

type_SpX_prop_tile <- neighbor_detail_summary |>
    ggplot(aes(x = Oligo.3_Astro, y = reference_SpX, fill = SpX_prop)) +
    geom_tile() +
    scale_fill_gradientn(colors = c(viridisLite::plasma(100))) + 
    theme_bw() +
    scale_y_discrete(limits=rev) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(type_SpX_prop_tile, filename = here(plot_dir, "Oligo3_nearest_astro_type_SpX_prop_tile.png"))



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

#### Neighbor Id vs. SpX ####

neighbor_SpX_prop <- neighbor_df_details |>
    count(BrNum, reference_SpX, neighbor_cell_type) |>
    group_by(BrNum, reference_SpX) |>
    mutate(prop = n/sum(n))


neighbor_SpX_prop_tile <- neighbor_SpX_prop |>
    ggplot(aes(x = neighbor_cell_type, y = reference_SpX, fill = prop)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("black", viridisLite::plasma(100))) +
    theme_bw() +
    scale_y_discrete(limits=rev)

ggsave(neighbor_SpX_prop_tile, filename = here(plot_dir, "Oligo3_nearest_astro_SpX_prop_tile.png"))

neighbor_SpX_n_tile <- neighbor_SpX_prop |>
    ggplot(aes(x = neighbor_cell_type, y = reference_SpX, fill = n)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("black", viridisLite::plasma(100))) +
    theme_bw() +
    scale_y_discrete(limits=rev)

ggsave(neighbor_SpX_n_tile, filename = here(plot_dir, "Oligo3_nearest_astro_SpX_n_tile.png"))

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
    ggplot(aes(x = distance, color = reference_SpX)) +
    geom_density() +
    scale_color_manual(values = metadata(spe)$SpX_colors) +
    geom_vline(xintercept = median_dist, color = "black", linetype = "dashed") +
    theme_bw()

ggsave(distance_density_SpX, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_density_SpX.png"))

distance_boxplot_SpX <- neighbor_df_details |>
    ggplot(aes(x = reference_SpX, y = distance, fill = reference_SpX)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = metadata(spe)$SpX_colors) +
    geom_hline(yintercept = median_dist, color = "gray25", linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(legend.position = "None")

ggsave(distance_boxplot_SpX, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_boxplot_SpX.png"))

distance_boxplot_donor <- neighbor_df_details |>
    ggplot(aes(x = BrNum, y = distance)) +
    geom_boxplot() +
    theme_bw() 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(distance_boxplot_donor, filename = here(plot_dir, "xenium_Oligo3_Astro_distance_boxplot_donor.png"))


#### Explore APOE expression ####

astro_apoe_expres <- plot_gene_express(spe[,spe$cell_type_broad == "Astro"], 
                                       genes = "APOE",
                                       category = "cell_type_anno",
                                       color_pal = cell_type_colors$anno)

ggsave(astro_apoe_expres, filename = here(plot_dir, "xenium_Astro_APOE_expres_violin.png"))

## Do nearest neighbor Astros express more APOE?

pd <- pd |>
    left_join(neighbor_df_details |>
                  group_by(neighbor_barcode) |>
                  slice_min(distance, with_ties = FALSE) |> 
                  select(cell_id = neighbor_barcode, reference_barcode, dist_class)) |>
    mutate(nn = ifelse(is.na(reference_barcode), "X", dist_class), # X = never a nearest neighbor
           cell_type_nn = paste0(cell_type_anno, "_", nn)) |>
    unique()

pd |> count(nn)

nrow(pd) == ncol(spe)

spe$cell_type_nn <- pd$cell_type_nn
spe$is_neighbor <- !is.na(pd$reference_barcode)

astro_apoe_expres <- plot_gene_express(spe[,spe$cell_type_broad == "Astro"], 
                                       genes = "APOE",
                                       category = "cell_type_nn")

ggsave(astro_apoe_expres, filename = here(plot_dir, "xenium_Astro_APOE_expres_violin_nn.png"))

astro_apoe_expres_nn_only <- plot_gene_express(spe[,spe$is_neighbor], 
                                       genes = "APOE",
                                       category = "cell_type_anno",
                                       color_pal = cell_type_colors$anno) +
    ggplot2::stat_summary( ## TODO make summary stat adjustable in plot_gene_express
        fun = median,
        geom = "crossbar",
        width = 0.3,
        color = "blue"
    )

ggsave(astro_apoe_expres_nn_only, filename = here(plot_dir, "xenium_Astro_APOE_expres_violin_nn_only.png"))


astro_apoe_expres_nn <- neighbor_df_details |>
    ggplot(aes(x=neighbor_cell_type, y=Astro_APOE, fill=neighbor_cell_type)) +
    geom_boxplot() +
    scale_fill_manual(values=cell_type_colors$anno) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(astro_apoe_expres_nn, filename = here(plot_dir, "xenium_Astro_APOE_expres_boxplot_nn.png"))

astro_apoe_expres_nn_dist <- neighbor_df_details |>
    mutate(neighbor_cell_type_dist = paste(neighbor_cell_type,  dist_class)) |>
    ggplot(aes(x=neighbor_cell_type_dist, y=Astro_APOE, fill=neighbor_cell_type, color=dist_class)) +
    geom_boxplot() +
    scale_fill_manual(values=cell_type_colors$anno) +
    scale_color_manual(values=c(near = "black", far = "grey70")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(astro_apoe_expres_nn_dist, filename = here(plot_dir, "xenium_Astro_APOE_expres_boxplot_nn_dist.png"))

#### APOE vs. Oligo distance #####

APOE_v_Oligo_dist <- neighbor_df_details |>
    ggplot(aes(distance, Astro_APOE)) +
    geom_point(size = 1, alpha = 0.1) +
    geom_smooth(method = "lm") +
    facet_wrap(~neighbor_cell_type) +
    theme_bw()

ggsave(APOE_v_Oligo_dist, filename = here(plot_dir, "xenium_Astro_APOE_v_Oligo_dist.png"))




# slurmjobs::job_single('20_xenium_Oligo3_Astro', create_shell = TRUE, memory = '20G', command = "Rscript 20_xenium_Oligo3_Astro.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

