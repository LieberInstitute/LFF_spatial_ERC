## Louise Huuki-Myers, May 2026
## Xenium visualizations 

#### Set Up ####

library("SpatialLIBD")
library("qs2")
library("here")
library("sessioninfo")
library("tidyverse")

data_dir <- here("processed-data", "21_Xenium", "15_xenium_vis")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "15_xenium_vis")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### load data ####
message(Sys.time(), " - Load SPE data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))


#### Layer markers ####

load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_SpD.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(gene %in% rownames(spe), MeanRatio >1) |>
    group_by(cellType.target) |>
    arrange(MeanRatio.rank) |>
    slice(1:5) |>
    select(gene, MeanRatio.rank, SpD = cellType.target) |>
    arrange(SpD)

marker_stats_top |> filter(grepl("WM", SpD))


map(c("ERMN"), function(gene){
    walk(unique(spe$BrNum), function(samp){
        vis_clus_class <- spatialLIBD::vis_gene(spe,
                                                sampleid = samp,
                                                geneid = gene,
                                                datatype = "Xenium",
                                                point_size = 1.5,
                                                guide_point_size = 3)
        
        ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_expres_%s_%s.png", gene, samp)), width = 12)
        
    })
})

    