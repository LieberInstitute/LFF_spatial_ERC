## Louise Huuki-Myers, May 2026
## Explore RCTD results

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
# library("BiocParallel")
library("tidyverse")
library("spatialLIBD")

# library("crumblr")
# library("variancePartition")

data_dir <- here("processed-data", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

## SET UP MULTICORE PARAM
# ncores <- 4  # match Sys.getenv('SLURM_CPUS_ON_NODE')
# bp <- MulticoreParam(workers = ncores) #or bp <- MulticoreParam(4)
# 

#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))


#### Add RCTD cell types to spe ####
message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

names(rctd_data@results)

head(rctd_data@results$results_df)

ncol(spe) == nrow(rctd_data@results$results_df)
identical(colnames(spe), rownames(rctd_data@results$results_df))

spe$cell_id <- colnames(spe)

colData(spe) <- cbind(colData(spe), rctd_data@results$results_df)

#### Explore RCTD results ####

table(is.na(spe$spot_class), spe$sum_gex < 100)

## from before sum_gex < 100 filter
#         FALSE   TRUE
# FALSE 450307      0
# TRUE       0  17209


table(rctd_data@results$results_df$spot_class)
# reject           singlet   doublet_certain doublet_uncertain 
# 5609            304401            122264             18033

table(rctd_data@results$results_df$first_type)

# Astro.1          Astro.2          Astro.3          Astro.4          Astro.5            Macro          Micro.1 
# 15125            18158              945             8954            42636             6285              291 
# Micro.2          Micro.3          Micro.4          Micro.5            OPC.1            OPC.2            OPC.3 
# 14780             2412              345            18178              120            11871               34 
# OPC.4            OPC.5          Oligo.1          Oligo.2          Oligo.3          Oligo.4          Oligo.5 
# 3814            11465            16888            11961            63536            37989             3599 
# Vasc.Endo          Vasc.PC        Vasc.VLMC         Excit.L2     Excit.L2_5.1     Excit.L2_5.2       Excit.L5.1 
# 36177            27006             8918            22124            10358             6601             8213 
# Excit.L5.2    Excit.L5_6_NP      Excit.L6_CT        Excit.L6b Inhib.Chandelier Inhib.Lamp5_Lhx6       Inhib.Pax6 
# 2341              983             5039             4974             1447             5284             1713 
# Inhib.Pvalb        Inhib.Sst        Inhib.Vip 
# 5247             4651             9845 

## spot class summary
rcdt_results_class_summary <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    count(sample_id, spot_class) |>
    group_by(sample_id) |>
    mutate(prop = n/sum(n))

## spot class summary plots
rctd_class_bar_plot <- rcdt_results_class_summary |>
    ggplot(aes(x = sample_id, y = n, fill = spot_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_bar_plot, filename = here(plot_dir, "xenium_rctd_class_bar_plot.png"))

rctd_class_prop_bar_plot <- rcdt_results_class_summary |>
    ggplot(aes(x = sample_id, y = prop, fill = spot_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_prop_bar_plot, filename = here(plot_dir, "xenium_rctd_class_prop_bar_plot.png"))

## cell type summary 
rcdt_results_singlets_summary <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    filter(spot_class == "singlet") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    group_by(sample_id, spot_class, cell_type_anno = first_type) |>
    summarise(xenium_n = n()) |>
    group_by(sample_id) |>
    mutate(xenium_singlet_prop = xenium_n/sum(xenium_n),
           type = "first_type")

rcdt_results_doublets <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    filter(spot_class == "doublet_certain") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    select(barcode, spot_class, sample_id, first_type, second_type) |>
    pivot_longer(!c(barcode, sample_id, spot_class), names_to = "type", values_to = "cell_type_anno") |> 
    bind_rows(rctd_data@results$results_df |> 
                  rownames_to_column("barcode") |>
                  filter(spot_class == "doublet_uncertain") |>
                  mutate(sample_id = gsub("_.*", "", barcode),
                         type = "first_type") |>
                  select(barcode, spot_class, sample_id, type,  cell_type_anno = first_type))

rcdt_results_doublets |> count(spot_class, type)

doublet_combos <- rcdt_results_doublets |> 
    filter(spot_class == "doublet_certain") |>
    group_by(barcode, sample_id) |> 
    summarise(doublet = paste0(sort(cell_type_anno), collapse = "-")) |>
    ungroup() |>
    count(doublet)

doublet_combos |> filter(grepl("Oligo.3", doublet)) |> arrange(-n)

# doublet                  n
# <chr>                <int>
# 1 Astro.5-Oligo.3       1714
# 2 Oligo.3-Excit.L6_CT   1615
# 3 Oligo.3-Excit.L6b     1608
# 4 Oligo.3-Excit.L2_5.1  1546
# 5 Oligo.3-Excit.L2_5.2  1022

rcdt_results_doublets_summary <- rcdt_results_doublets |>    
    group_by(sample_id, spot_class, type, cell_type_anno) |>
    summarise(xenium_n = n())  |>
    group_by(sample_id) |>
    mutate(xenium_doublet_prop = xenium_n/sum(xenium_n))

rcdt_results_summary <- rcdt_results_singlets_summary |> 
    bind_rows(rcdt_results_doublets_summary)  |>
    group_by(sample_id, cell_type_anno) |>
    summarise(xenium_n = sum(xenium_n)) |>
    group_by(sample_id) |>
    mutate(xenium_prop = xenium_n/sum(xenium_n))


cell_type_class_summary <- rcdt_results_singlets_summary |> 
    bind_rows(rcdt_results_doublets_summary) |>
    group_by(cell_type_anno, spot_class, type) |>
    summarise(xenium_n = sum(xenium_n)) |>
    group_by(cell_type_anno) |> 
    mutate(prop = xenium_n/sum(xenium_n),
           rctd_class = ifelse(spot_class == "doublet_certain", paste0("doublet_c_", gsub("_type", "", type)) , as.character(spot_class)))

rctd_class_ct_bar_plot <- cell_type_class_summary |>
    ggplot(aes(x = cell_type_anno, y = xenium_n, fill = rctd_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_ct_bar_plot, filename = here(plot_dir, "xenium_rctd_class_ct_bar_plot.png"), width = 9)

rctd_class_ct_prop_bar_plot <- cell_type_class_summary |>
    ggplot(aes(x = cell_type_anno, y = prop, fill = rctd_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_ct_prop_bar_plot, filename = here(plot_dir, "xenium_rctd_class_ct_prop_bar_plot.png"), width = 9)

write_csv(cell_type_class_summary, file = here(data_dir, "RCTD_cell_type_class_summary.csv"))

#Save SPE with additional data
message(Sys.time(), " - Saving cleaned SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_cell_types.qs2"))

# spe <- qs_read(here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2"))

source(here("code", "utils", "my_plot_reduced_dim.R"))

#### plot ####
## categorical
walk2(c("first_type"), list(cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat", color_pal = .y))
walk2(c("first_type"), list(cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", color_pal = .y))

walk2(c("second_type"), list(cell_type_colors$anno), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", color_pal = .y))

walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat"))
walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat", facet = TRUE))

walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat"))
walk(c("spot_class"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat", facet = TRUE))


#### doublets ####

rctd_tb_doublet_counts <- rctd_tb |>
    filter(spot_class == "doublet_certain") |>
    count(first_type, second_type) |>
    as_tibble()

rctd_tb_doublet_counts |> arrange(-n) 

rctd_tb_doublet_counts |> filter(first_type == "Oligo.3") |> arrange(-n) 
    

rctd_tb_doublet_tile <- rctd_tb_doublet_counts |>
    ggplot(aes(first_type, second_type, fill = n)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(rctd_tb_doublet_tile, filename = here(plot_dir, "xenium_rctd_doublet_tile.png"))

#### Spatial plots ####

vis_clus_class <- vis_clus(spe, 
                           sampleid = "Br1556",
                           clustervar = "spot_class",
                           datatype = "Xenium", 
                           point_size = 1.5,
                           # alpha = 0.5,
                           colors = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
                           guide_point_size = 2)

ggsave(vis_clus_class, filename = here(plot_dir, "Xenium_vis_spot_class_Br1556.png"), width = 12)


vis_clus_first_type <- vis_clus(spe, 
                           sampleid = "Br1556",
                           clustervar = "first_type",
                           datatype = "Xenium", 
                           point_size = 1.5,
                           colors = cell_type_colors$anno,
                           guide_point_size = 2)

ggsave(vis_clus_first_type, filename = here(plot_dir, "Xenium_vis_first_type_Br1556.png"), width = 12)

vis_clus_first_type_singlet <- vis_clus(spe[, spe$spot_class == "singlet"], 
                                sampleid = "Br5415",
                                clustervar = "first_type",
                                datatype = "Xenium", 
                                point_size = 1.5,
                                colors = cell_type_colors$anno,
                                guide_point_size = 2)

ggsave(vis_clus_first_type_singlet, filename = here(plot_dir, "Xenium_vis_first_type_singlet_Br5415.png"), width = 12)

vis_clus_first_type_Oligo.3 <- vis_clus(spe[, spe$spot_class == "singlet" & spe$first_type == "Oligo.3"], 
                                        sampleid = "Br1556",
                                        clustervar = "first_type",
                                        datatype = "Xenium", 
                                        point_size = 1.5,
                                        colors = cell_type_colors$anno,
                                        guide_point_size = 2)

ggsave(vis_clus_first_type_Oligo.3, filename = here(plot_dir, "Xenium_vis_first_type_Oligo.3_Br1556.png"), width = 12)

vis_clus_first_type_Oligo.1 <- vis_clus(spe[, spe$spot_class == "singlet" & spe$first_type == "Oligo.1"], 
                                        sampleid = "Br1556",
                                        clustervar = "first_type",
                                        datatype = "Xenium", 
                                        point_size = 1.5,
                                        colors = cell_type_colors$anno,
                                        guide_point_size = 2)

ggsave(vis_clus_first_type_Oligo.1, filename = here(plot_dir, "Xenium_vis_first_type_Oligo.1_Br1556.png"), width = 12)




#### SC dot plots ####


# slurmjobs::job_single('10_xenium_cell_types', create_shell = TRUE, memory = '100G', command = "Rscript 10_xenium_cell_types.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

