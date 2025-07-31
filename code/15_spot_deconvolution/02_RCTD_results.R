## Louise Huuki-Myers, July 2025
## Process and plot RCTD spot deconvolution data

#### set up ####
library("ggplot2")
library("SpatialExperiment")
library("spacexr")
library("HDF5Array")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("tidyverse")

plot_dir <- here("plots", "15_spot_deconvolution", "02_RCTD_results")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "15_spot_deconvolution", "02_RCTD_results")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)

cell_type_broad_levels <- names(cell_type_colors$broad)
cell_type_broad_levels <- cell_type_broad_levels[cell_type_broad_levels != "Other"]

#### load results ####
results_spe <- readRDS(here("processed-data", "15_spot_deconvolution", "01_runRCTD", "RCDT_est_prop.rds"))
colData(results_spe)

table(results_spe$sample_id)

# Assays have cell types as rows and pixels as columns

assay(results_spe, "weights")

summary(colSums(assay(results_spe, "weights")))

sort(rowSums(assay(results_spe, "weights")))

## Load HD5F spe for colData
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
message("spots failed RCDT: ", ncol(spe) - ncol(results_spe)) #40

spe <- spe[,colnames(results_spe)]

sample_info <- colData(spe)[,c("BrNum","SpD", "APOE", "CNmask_dark_blue")]

sample_ids(results_spe)

identical(colnames(spe), colnames(results_spe))

colData(results_spe) <- colData(results_spe)
colData(results_spe) <- cbind(sample_info, colData(results_spe))

#### classification ####
classification_df <- data.frame(
    pixel = colnames(assay(results_spe)),
    spot_class = colData(results_spe)$spot_class,
    first_type = colData(results_spe)$first_type,
    second_type = colData(results_spe)$second_type
)

#### plot composition for each spot ####

results_spe[,results_spe$BrNum == "Br5212"]

pdf(here(plot_dir, "test.pdf"))
plotAllWeights(results_spe[,results_spe$BrNum == "Br5212"], 
               assay_name = "weights",
               r = 100, 
               lwd = 0
               )
dev.off()

pdf(here(plot_dir, "test2.pdf"))
plotCellTypeWeight(
    results_spe[,results_spe$BrNum == "Br5212"], 
    "Oligo.3", 
    "weights_full",
    size = 2, 
    stroke = 0.2, 
    title = "\nDensity of Oligo.3\n"
)
dev.off()



#### replace our assays with RCDT results ####


results_spe2 <- SpatialExperiment(sample_id = spe$sample_id,
                  colData = colData(spe),
                  assays = list(weights = assay(results_spe, "weights"),
                                weights_unconfident = assay(results_spe, "weights_unconfident"),
                                weights_full = assay(results_spe, "weights_full")),
                  spatialCoords = spatialCoords(spe),
                  imgData = imgData(spe)
                  )

assay(results_spe2, "cell_counts") <- ceiling(assay(results_spe2, "weights") * results_spe2$CNmask_dark_blue)
assay(results_spe2, "cell_counts")[1:5, 1:5]

test_vis_gene <- vis_gene(results_spe2, 
         sampleid = "Br5212",
         # geneid = c("Oligo.1","Oligo.2", "Oligo.3","Oligo.4","Oligo.5"),
         geneid = "Astro.5",
         assayname = "weights",
         cont_colors = viridisLite::rocket(10, direction = -1))

ggsave(test_vis_gene, filename = here(plot_dir, "test_vis_gene.png"))


#### long data ####

sample_info <- colData(spe)[,c("key", "BrNum","SpD", "APOE", "CNmask_dark_blue")] |>
    as.data.frame() 

rcdt_long <- assay(results_spe2, "weights") |>
    as.matrix() |>
    reshape2::melt() |>
    dplyr::rename(cell_type_fine = Var1, key = Var2, weights = value) |>
    as_tibble() |>
    left_join(sample_info) |>
    mutate(cell_type_broad = factor(jaffelab::ss(as.character(cell_type_fine), "\\."),
                                    levels = cell_type_broad_levels))

rcdt_long |> count(cell_type_broad)

rcdt_summary <- rcdt_long |>
    group_by(cell_type_broad, cell_type_fine, BrNum, SpD) |>
    summarise(max = max(weights),
              mean = mean(weights),
              median = median(weights),
              min = min(weights))


rcdt_violin <- rcdt_long |>
    filter(grepl("Excit", cell_type_fine))|>
    ggplot(aes(x = SpD, y = weights, color = SpD)) +
    geom_violin() +
    facet_wrap(~cell_type_fine, ncol = 1) +
    theme_bw() +
    scale_color_manual(values = SpD_colors) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(rcdt_violin, filename = here(plot_dir, "rcdt_violin.png"), height = 10)


pdf(here(plot_dir, "rcdt_mean_boxplot.pdf"))
map(cell_type_broad_levels, ~rcdt_summary |>
        filter(cell_type_broad == .x)|>
        ggplot(aes(x = SpD, y = mean, color = SpD)) +
        geom_boxplot() +
        facet_wrap(~cell_type_fine, ncol = 1) +
        theme_bw() +
        scale_color_manual(values = SpD_colors) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    )
dev.off()

ggsave(rcdt_mean_boxplot, filename = here(plot_dir, "rcdt_mean_boxplot.png"), height = 10)




