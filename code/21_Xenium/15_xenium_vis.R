## Louise Huuki-Myers, May 2026
## Xenium visualizations 

#### Set Up ####

library("spatialLIBD")
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

## local
# spe <- qs2::qs_read(here::here("code","23_spatialLIBD_app_Xenium","spe_xenium_app.qs2"))

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


map(c("MBP", "MBP", "ERMN"), function(gene){
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


SpX_markers <- list(
    `Inhib~Sp9X5`  = "VIP",        # canonical inhibitory marker, highest logFC (2.39) of specific markers
    `L1~Sp9X6`     = "SFRP2",      # highest logFC (4.02) of any L1 gene, astrocyte/L1 specific
    `L1~Sp9X7`     = "FGFR3",      # clean homeostatic astrocyte marker, strong stat (8.01)
    `L2.3~Sp9X4`   = "ADAMTS3",    # top stat (13.44), validated L2 specific marker in your panel
    `L5~Sp9X1`     = "PCSK1",      # L5_ET base probe, strong logFC (2.08), projection neuron specific
    `L6~Sp9X9`     = "RXFP1",      # top 3 stat (13.11), highest logFC (2.94), clean L6 specificity
    `LD~Sp9X8`     = "CNDP1",      # top stat for LD, Oligo base probe
    `Vasc~Sp9X3`   = "PECAM1",     # canonical endothelial marker
    `WM~Sp9X2`     = "ERBB3"       # Oligo base probe, strong stat (11.31), clean WM biology
)

#### Vis clus ####
SpX_colors <- metadata(spe)$SpX_colors

samp = "Br6538"
vis_clus_class <- spatialLIBD::vis_clus(spe,
                                        sampleid = samp,
                                        clustervar = "SpX",
                                        datatype = "Xenium",
                                        point_size = 1.5,
                                        # alpha = 0.5,
                                        colors = SpX_colors,
                                        guide_point_size = 3)

ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_SpX_%s.png", samp)), width = 12, height = 11)


vis_grid_clus(spe,
              clustervar = "SpX",
              pdf_file = here(plot_dir, "xenium_orient_test_SpX.pdf"),
              colors = SpX_colors,
              point_size = 1,
              sample_order = sort(unique(spe$sample_id)),
              datatype = "Xenium")

    