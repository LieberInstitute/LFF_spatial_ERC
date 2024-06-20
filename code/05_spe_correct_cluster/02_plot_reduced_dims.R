library("SpatialExperiment")
library("here")
library("sessioninfo")
library("scran")
library("scater")
library("Harmony")

plot_dir <- here("plots", "05_spe_correct_cluster", "02_plot_reduced_dims")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "02_build_spe", "01_preprocess_Harmony")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

spe <- readRDS(here("processed-data", "02_build_spe", "spe.rds"))

#### plot Uncorrected Data ####
qc_colors <- c(fold = "#4BA402", out_edge = "#097FE0", vessel ="#BB1EF4", None = "#CCCCCC40")
scran_colors <- c(`TRUE` = "red", `FALSE` = "#CCCCCC40")


##umap
umap_sample <- ggcells(spe, mapping = aes(x = UMAP.1, y = UMAP.2, color = sample_id))+
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

ggsave(umap_sample, filename = here(plot_dir, "UMAP_sample_id.png"))

umap_sample_facet <- ggcells(spe, mapping = aes(x = UMAP.1, y = UMAP.2, color = sample_id))+
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    facet_wrap(~sample_id) +
    theme(legend.position = "None")

ggsave(umap_sample_facet, filename = here(plot_dir, "UMAP_sample_id_facet.png"))

umap_qc_anno <- ggcells(spe, mapping = aes(x = UMAP.1, y = UMAP.2, color = qc_anno))+
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    scale_color_manual(values = qc_colors) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

ggsave(umap_qc_anno , filename = here(plot_dir, "UMAP_qc_anno.png"))

umap_scran <- ggcells(spe, mapping = aes(x = UMAP.1, y = UMAP.2, color = scran_discard))+
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    scale_color_manual(values = scran_colors) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

ggsave(umap_scran , filename = here(plot_dir, "UMAP_scran.png"))

umap_sum_umi <- ggcells(spe, mapping = aes(x = UMAP.1, y = UMAP.2, color = sum_umi))+
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    scale_color_viridis() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

ggsave(umap_scran , filename = here(plot_dir, "umap_sum_umi.png"))
