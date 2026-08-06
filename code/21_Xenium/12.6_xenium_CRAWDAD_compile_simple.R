## Louise Huuki-Myers, May 2026
## Compile & plot CRAWDAD results from all samples
## Adapted from https://github.com/LieberInstitute/Habenula_Visium/blob/01fc3561339376fb5333163fbbcbb6f16d12ddaf/code/09_HD_cell_level/16_crawdad_plot.R

#### Set Up ####

library("here")
library("sessioninfo")
library("crawdad")
library("tidyverse")
library("spatialLIBD")
library("scales")
library("qs2")
library("patchwork")

data_dir <- here("processed-data", "21_Xenium", "12.6_xenium_CRAWDAD_compile_simple")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "12.6_xenium_CRAWDAD_compile_simple")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)
cell_type_levels <- names(cell_type_colors$anno)

#### get layer organization ####
load(here("processed-data", "SpX_colors.Rdata"))

## cell type SpX organization
cell_v_SpX_prop_long <- read.csv(file = here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "cell_v_SpX_prop_long.csv"), row.names = 1)

cell_v_SpX_max <- cell_v_SpX_prop_long |> 
    mutate(SpX = factor(SpX, levels = names(SpX_colors))) |>
    group_by(cell_type_anno) |>
    slice_max(n_cell) |> 
    arrange(SpX, cell_type_anno)

cell_v_SpX_max |> dplyr::count(SpX)

#### Get CRAWDAD data ####

crawdad_data_dir <- here("processed-data", "21_Xenium", "12_xenium_CRAWDAD")

crawdad_data <- map_dfr(list.files(crawdad_data_dir, full.names = TRUE), ~read_csv(.x, show_col_types = FALSE) |>
                            mutate(file = basename(.x)))  |>
    mutate(
        id = BrNum,
        reference = factor(reference, levels = cell_type_levels),
        neighbor = factor(neighbor, levels = cell_type_levels)
    ) |> 
    filter(region == "ALL", scale <=50)

## Missing data - resubmit jobs
# crawdad_data |> count(region)
# 
# missing_brnum <- crawdad_data |> filter(is.na(region)) |> count(BrNum) |> pull(BrNum)
# 
# logs_tb <- tibble(log_fn = list.files(here("code", "21_Xenium", "logs"), pattern = "12_xenium_CRAWDAD")) |>
#     separate(log_fn, into = c(NA, NA, NA, "BrNum", "region", "jobid"), remove = FALSE)
# 
# missing_brnum_region <- tidyr::expand_grid(BrNum = missing_brnum, region = c("ALL", "Vasc", "GM", "WM")) |>
#     anti_join(crawdad_data |> filter(BrNum %in% missing_brnum) |> count(BrNum, region)) |>
#     left_join(logs_tb) |>
#     mutate(jobid = as.numeric(jobid))|>
#     arrange(jobid)
# 
# cat(missing_brnum_region$jobid, sep = ", ")

any(is.na(crawdad_data$region))

crawdad_data |>
    dplyr::count(is.na(Z), 
                 reference) |>
    arrange(-n)


#   Get the Z-score significance threshold (same in all samples)
(z_sig <- crawdad_data |> filter(BrNum == "Br1039") |> correctZBonferroni())
# [1] 4.12

crawdad_data |>
    filter(reference == 'Oligo.3', neighbor =='Oligo.1') |>
    group_by(BrNum, region, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(BrNum, region, neighbor, reference) |>
    filter(scale == min(scale)) |>
    arrange(Z)


trend_test <- crawdad_data |> 
    filter(reference == 'Astro.4', grepl("Vasc", neighbor)) |>
    group_by(id = BrNum, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |> 
    vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig) +
    facet_wrap(~neighbor)

ggsave(trend_test, filename  = here(plot_dir, "xenium_CRAWDAD_trend_Astro.4_Vasc.png"))

trend_test <- crawdad_data |> 
    filter(reference == 'Oligo.5', grepl("Vasc", neighbor)) |>
    group_by(id = BrNum, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |> 
    vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig) +
    facet_wrap(~neighbor)

ggsave(trend_test, filename  = here(plot_dir, "xenium_CRAWDAD_trend_Oligo.5_Vasc.png"))


trend_test <- crawdad_data |> 
    filter(neighbor == 'Oligo.3', grepl("Astro", reference)) |>
    group_by(id = BrNum, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |> 
    vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig) +
    facet_wrap(~reference)

ggsave(trend_test, filename  = here(plot_dir, "xenium_CRAWDAD_trend_Oligo.3_Astro.png"))

crawdad_data |> dplyr::count(is.na(Z))
crawdad_data |> filter(is.na(Z)) |> dplyr::count(BrNum) |> arrange(-n)

#### Summarize across samples and regions ####
min_num_signif = 4

crawdad_data_summary <- crawdad_data |>
    filter(!is.na(Z)) |> ## filter missing scores
    group_by(region, BrNum, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(region, BrNum, neighbor, reference) |>
    filter(scale == min(scale)) |>
    #   Retain pairs where all samples all signs of Z scores agree across
    #   samples, and significance is achieved in some sufficient number of
    #   samples 
    group_by(region, reference, neighbor) |>
    filter(mean(Z > 0) >= .8 | mean(Z < 0) >= 0.8) |> ## allows 20% disagreement 
    filter(n() >= min_num_signif) |>
    #   Take the mean Z-score and scale across samples
    group_by(region, neighbor, reference) |>
    summarize(scale = mean(scale), Z = mean(Z)) |>
    ungroup() |>
    #   Cap Z-score at twice the magnitude of the significance threshold 
    mutate(Z_real = Z,
           Z = sign(Z) * pmin(abs(Z), z_sig * 2))

write_csv(crawdad_data_summary, file = here(data_dir, "ERC_xenium_crawdad_data_summary_All_max50.csv"))

## check pairs of interest
crawdad_data_summary |>
    # filter(reference == "Oligo.3", neighbor == "Excit.L2")
    # filter(reference == "Oligo.3", grepl("Astro", neighbor))
    filter(neighbor == "Oligo.3", grepl("Astro", reference))

crawdad_data_summary |>
    filter(reference == "Oligo.3") |>
    arrange(-Z_real) 

crawdad_data_summary |>
    filter(reference == "Astro.4") |>
    arrange(-Z_real) |>
    print(n = 32)

crawdad_data_summary |>
    filter(reference == "Oligo.3") |>
    # filter(neighbor == "Oligo.3") |> 
    arrange(-Z_real) |>
    print(n=30)

crawdad_data_summary |>
    # filter(grepl("Vasc", reference) & grepl("Astro", neighbor) ) |> 
    filter(grepl("Astro", reference) & grepl("Vasc", neighbor) ) |> 
    arrange(-Z_real) |>
    print(n=30)

crawdad_data_summary |>
    filter(neighbor == "Astro.4") |> 
    arrange(-Z_real) |>
    print(n=30)

crawdad_data_summary |>
    filter(reference == "Oligo.5")

crawdad_data_summary |>
    # filter(reference == "Astro.4") |>
    filter(neighbor == "Astro.4") |>
    arrange(-Z_real) 


crawdad_data_summary |>
    filter(reference == "Astro.1") |>
    group_by(reference) |>
    slice_max(abs(Z_real))

crawdad_data_summary_top <- crawdad_data_summary |>
    group_by(reference) |>
    slice_max(abs(Z_real)) |>
    slice_min(scale)

pdf(here(plot_dir, "xenium_CRAWDAD_trend_plots.pdf"))

crawdad_data_summary_top |>
    dplyr::select(ref = reference, nei = neighbor, summary_scale = scale, Z_real)|>
    pmap(function(ref, nei, summary_scale, Z_real){
        
        trend_plot <- crawdad_data |> 
            filter(reference == ref, neighbor == nei)  |>
            group_by(id = BrNum, neighbor, scale, reference) |>
            summarize(Z = mean(Z)) |>
            vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig) +
            labs(title = sprintf("Z=%.2f, scale=%.2f", Z_real, summary_scale))
        
        return(trend_plot)
        
    })

dev.off()

pdf(here(plot_dir, "xenium_CRAWDAD_trend_plots_Oligo.3.pdf"))

crawdad_data_summary |>
    filter(reference == "Oligo.3") |>
    dplyr::select(ref = reference, nei = neighbor, summary_scale = scale, Z_real)|>
    pmap(function(ref, nei, summary_scale, Z_real){
        
        trend_plot <- crawdad_data |> 
            filter(reference == ref, neighbor == nei)  |>
            group_by(id = BrNum, neighbor, scale, reference) |>
            summarize(Z = mean(Z)) |>
            vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig) +
            labs(title = sprintf("Z=%.2f, scale=%.2f", Z_real, summary_scale))
        
        return(trend_plot)
        
    })

dev.off()

pdf(here(plot_dir, "xenium_CRAWDAD_trend_plots_Astro.4.pdf"))

crawdad_data_summary |>
    filter(reference == "Astro.4") |>
    dplyr::select(ref = reference, nei = neighbor, summary_scale = scale, Z_real)|>
    pmap(function(ref, nei, summary_scale, Z_real){
        
        trend_plot <- crawdad_data |> 
            filter(reference == ref, neighbor == nei)  |>
            group_by(id = BrNum, neighbor, scale, reference) |>
            summarize(Z = mean(Z, na.rm = TRUE)) |>
            vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = z_sig) +
            labs(title = sprintf("Z=%.2f, scale=%.2f", Z_real, summary_scale))
        
        return(trend_plot)
        
    })

dev.off()

## TODO check this
crawdad_data |>
    filter(reference == "Astro.4", 
           grepl("Vasc", neighbor))

#### custom dot plot ####

# pivot to matrix for clustering
z_mat <- crawdad_data |>
    group_by(reference, neighbor) |>
    summarise(Z=mean(Z), .groups="drop") |>
    pivot_wider(names_from=reference, values_from=Z, values_fill=0) |>
    column_to_rownames("neighbor") |>
    as.matrix()

# cluster both axes
neighbor_order <- hclust(dist(z_mat))$labels[hclust(dist(z_mat))$order]
reference_order <- hclust(dist(t(z_mat)))$labels[hclust(dist(t(z_mat)))$order]


result_df <- crawdad_data |> 
    left_join(cell_v_SpX_max |> 
                  select(neighbor_SpX = SpX, 
                         neighbor = cell_type_anno)
              )|> 
    left_join(cell_v_SpX_max |> 
                  select(reference_SpX = SpX, 
                         reference = cell_type_anno)
              ) |>
    mutate(
        neighbor=fct_relevel(neighbor, neighbor_order),
        reference=fct_relevel(reference, neighbor_order) ## for now - ref would not cluster due to NA
    ) 

custom_dotplot = function(result_df, z_sig) {
    p = result_df |>
        ggplot(
        aes(x = reference, y = neighbor, color = Z)
    ) +
        geom_point() +
        scale_color_gradient2(low="blue", mid="#CECECE", high="red", midpoint=0) +
        # scale_radius(
        #     trans = 'reverse',
        #     breaks = seq(
        #         min(result_df$scale), max(result_df$scale), length.out = 3
        #     ),
        #     range = c(2, 8)
        # ) +
        # facet_grid(neighbor_SpX~reference_SpX, scales = "free") +
        # coord_fixed() +
        # theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    
    return(p)
}


walk(c("ALL", "Vasc", "GM", "WM"), function(reg){
    dotplot <- custom_dotplot(crawdad_data_summary, z_sig, reg = reg)
        ggsave(dotplot, filename = here(plot_dir, sprintf('erc_xenium_CRAWDAD_dot_plot-%s.png', reg)), height = 12 , width = 12)
    
})


crawdad_data_summary <- crawdad_data_summary |> 
    mutate(referencer = factor(reference, levels = cell_v_SpX_max$cell_type_anno),
           neighbor = factor(neighbor, levels = cell_v_SpX_max$cell_type_anno))

walk(c("ALL", "Vasc", "GM", "WM"), function(reg){
    dotplot <- custom_dotplot(crawdad_data_summary, z_sig, reg = reg)
    ggsave(dotplot, filename = here(plot_dir, sprintf('erc_xenium_CRAWDAD_dot_plot_layer-%s.png', reg)), height = 12 , width = 12)
    
})


library(ComplexHeatmap)

Heatmap(z_mat, cluster_columns = FALSE)


#### Visualize in Sections ####

## load spe  data 
message(Sys.time(), " - Load SPE data")
spe <- qs2::qs_read(here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2"))
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

spe <- spe[,spe$spot_class == "singlet"]


region_tb <- list("WM"= c('WMtz~SpX8', 'WM~SpX2'),
     "GM"= c("L2.3~SpX4", "Inhib~SpX5", "L5~SpX1", "L6~SpX9"),
     "Vasc"= c("Vasc~SpX3", "L1~SpX6", "L1~SpX7"))  |>
    enframe(name="group", value="cluster") |>
    unnest(cluster) |>
    dplyr::rename(region = group, SpX = cluster)


spe$region <- region_tb$region[match(spe$SpX, region_tb$SpX)]

table(spe$SpX, spe$region)

spe_test <- spe[,spe$BrNum == "Br1039"]

summary(spatialCoords(spe_test)[,"x_centroid"])
summary(spatialCoords(spe_test)[,"y_centroid"])


spe_test[, spatialCoords(spe_test)[,"y_centroid"] < 850]

vis_clus_region <- vis_clus(spe,
         sampleid = "Br1039",
        clustervar = "region",
        datatype = "Xenium",
        # colors = cell_type_colors$anno, 
        point_size = 1) +
    theme_bw() +
    geom_hline(yintercept = min(spatialCoords(spe_test)[,"y_centroid"]) + c(50, 100, 200, 300, 400, 500, 1000, 5000))

ggsave(vis_clus_region, filename = here(plot_dir, "vis_clus_region_scales.png"), height = 10, width = 10)


plot_ref_neighbor_pair <- function(spe, sampleid = "Br1039", ref, neighbor){
    
    spe$ref_neighbor <- spe$cell_type_anno
    
    spe$ref_neighbor[!spe$ref_neighbor %in% c(ref, neighbor)] <- NA
    spe$ref_neighbor <- droplevels(spe$ref_neighbor)
    
    ref_neighbor_colors <- c("red", "blue")
    names(ref_neighbor_colors) <- c(ref, neighbor)
    
    cell_type_vis_clus <- vis_clus(spe,
                                   sampleid = sampleid,
                                   clustervar = "ref_neighbor",
                                   datatype = "Xenium",
                                   colors = ref_neighbor_colors, 
                                   point_size = 1.5)
    
    ggsave(cell_type_vis_clus, filename = here(plot_dir, sprintf("vis_clus_%s_%s_%s.png", sampleid, ref, neighbor)), height = 10, width = 10)
    
}

crawdad_data_summary |>
    filter(reference == "Oligo.3") |>
    arrange(Z_real)

plot_ref_neighbor_pair(spe, ref = "Astro.1", neighbor = "Vasc.Endo")
plot_ref_neighbor_pair(spe, ref = "Astro.4", neighbor = "Vasc.Endo")
plot_ref_neighbor_pair(spe, ref = "Oligo.3", neighbor = "Oligo.1")

plot_ref_neighbor_pair(spe, ref = "Oligo.3", neighbor = "Astro.1")
plot_ref_neighbor_pair(spe, sampleid = "Br5941", ref = "Oligo.3", neighbor = "Astro.1")
plot_ref_neighbor_pair(spe, sampleid = "Br6098", ref = "Oligo.3", neighbor = "Astro.1")


plot_ref_neighbor_pair(spe, sampleid = "Br5460", ref = "Oligo.3", neighbor = "Astro.1")
plot_ref_neighbor_pair(spe, sampleid = "Br5460", ref = "Oligo.3", neighbor = "Astro.2")
plot_ref_neighbor_pair(spe, sampleid = "Br5460", ref = "Oligo.3", neighbor = "Astro.3")

plot_ref_neighbor_pair(spe, sampleid = "Br5460", ref = "Oligo.3", neighbor = "Vasc.Endo")
plot_ref_neighbor_pair(spe, sampleid = "Br5460", ref = "Oligo.5", neighbor = "Vasc.Endo")
plot_ref_neighbor_pair(spe, ref = "Oligo.3", neighbor = "Oligo.4")


plot_ref_neighbor_pair(spe, ref = "Oligo.1", neighbor = "Excit.L5.2")
plot_ref_neighbor_pair(spe, ref = "Oligo.3", neighbor = "Excit.L5.2")
