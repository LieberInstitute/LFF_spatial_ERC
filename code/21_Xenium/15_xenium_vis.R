## Louise Huuki-Myers, May 2026
## Xenium visualizations 

#### Set Up ####

library("spatialLIBD")
library("qs2")
library("here")
library("sessioninfo")
library("tidyverse")
library("patchwork")

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

load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_vSpD.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(gene %in% rownames(spe), MeanRatio >1) |>
    group_by(cellType.target) |>
    arrange(MeanRatio.rank) |>
    slice(1:5) |>
    select(gene, MeanRatio.rank, vSpD = cellType.target) |>
    arrange(vSpD)

marker_stats_top |> filter(grepl("WM", vSpD))


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

#### Vis clus ####
SpX_colors <- metadata(spe)$SpX_colors

samp = "Br6538"
vis_clus_class <- spatialLIBD::vis_clus(spe,
                                        sampleid = samp,
                                        clustervar = "xSpD",
                                        datatype = "Xenium",
                                        point_size = 1.5,
                                        # alpha = 0.5,
                                        colors = SpX_colors,
                                        guide_point_size = 3)

ggsave(vis_clus_class, filename = here(plot_dir, sprintf("Xenium_SpX_%s.png", samp)), width = 12, height = 11)


vis_grid_clus(spe,
              clustervar = "xSpD",
              pdf_file = here(plot_dir, "xenium_orient_test_SpX.pdf"),
              colors = SpX_colors,
              point_size = 1,
              sample_order = sort(unique(spe$sample_id)),
              datatype = "Xenium")

#### Build Fig 1 vis_clus grid ####

# representative sections:
# Br1556	AA	E2
# Br2582	EA	E2 (other sample Br6538 is the aligment issue)
# Br5460	AA	E4 (Br1039 also okay)
# Br5517	EA	E4	(adjacnet gyrus changes shapes between secitons)

rep_sections <- c("Br1556", "Br2582", "Br5460", "Br5517")

#### Spot plots for representative sections ####

single_vis_clus <- vis_clus(
    spe = spe,
    point_size = 1.7,
    colors = metadata(spe)$SpX_colors,
    sampleid = "Br5517",
    clustervar = "xSpD",
    spatial = FALSE,
    guide_point_size = 5,
    datatype = "Xenium"
)

ggsave(single_vis_clus, filename = here(plot_dir, "xenium_vis_clus_Br5517.png")) 


xenium_row_plots <- map(rep_sections, function(s) {
    vis_clus_plot <- vis_clus(
        spe = spe,
        point_size = 1.2,
        colors = metadata(spe)$SpX_colors,
        sampleid = s,
        clustervar = "xSpD",
        datatype = "Xenium",
        guide_point_size = 5
    ) +
        labs(title = NULL)  +
        theme(
            legend.position = "None", ## using heat maps label colors
            axis.title.x = element_blank(),
            text = element_text(size = 14)
        )
    
    return(vis_clus_plot)
})

## add legend to last plot
xenium_row_plots[[4]] <- xenium_row_plots[[4]] + theme(legend.position = "right")

xenium_row <- Reduce("|", xenium_row_plots)
# ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.pdf"), width = 18, height = 9)

ggsave(xenium_row, filename = here(plot_dir, "vis_SpX_rep_sections.png"), width = 18, height = 4.5)

#### plot vis sections for Visium data ####

message(Sys.time(), " - Load HDF5 SPE")
spe_vis <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))

visium_row_plots <- map(rep_sections, function(s) {
    vis_clus_plot <- vis_clus(
        spe = spe_vis,
        point_size = 1.3,
        colors = metadata(spe_vis)$SpD_colors,
        sampleid = s,
        clustervar = "vSpD",
        guide_point_size = 5
    ) +
        labs(title = s)  +
        theme(
            legend.position = "None", ## using heat maps label colors
            axis.title.x = element_blank(),
            text = element_text(size = 14)
        )
    
    return(vis_clus_plot)
}) 

## add legend to last plot
visium_row_plots[[4]] <- visium_row_plots[[4]] + theme(legend.position = "right")

visium_row <- Reduce("|", visium_row_plots)
# ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.pdf"), width = 18, height = 9)

ggsave(visium_row, filename = here(plot_dir, "vis_SpD_rep_sections.png"), width = 18, height = 4.5)


xen_vis_rep_grid <- visium_row / xenium_row
ggsave(xen_vis_rep_grid, filename = here(plot_dir, "vis_SpD_SpX_rep_sections.png"), width = 18, height = 9)

## plot APOE in both Xenium and Visium
rownames(spe) <- rowData(spe)$gene_name
APOE_v_xSpD <- DeconvoBuddies::plot_gene_express(spe_pb, genes = "APOE", category = "xSpD", color_pal = metadata(spe)$SpX_colors)

ggsave(APOE_v_xSpD, filename = here(plot_dir, "ERC_Xenium_SpD_express_violin_APOE.png"), height = 4, width = 6)


rownames(spe_vis) <- rowData(spe_vis)$gene_name
APOE_v_vSpD <- DeconvoBuddies::plot_gene_express(spe_vis, genes = "APOE", category = "vSpD", color_pal = metadata(spe_vis)$SpD_colors)


ggsave(APOE_v_vSpD + APOE_v_vSpD, filename = here(plot_dir, "ERC_both_SpD_express_violin_APOE.png"), height = 4, width = 6)

# slurmjobs::job_single('15_xenium_vis', create_shell = TRUE, memory = '50G', command = "Rscript 15_xenium_vis.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()