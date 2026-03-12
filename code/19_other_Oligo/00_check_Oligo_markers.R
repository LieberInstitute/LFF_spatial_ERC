## Louise Huuki-Myers, Match 2025
## Check Oligo marker genes (data drven or lit)

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("DeconvoBuddies")

data_dir <- here("processed-data", "19_other_Oligo", "00_check_Oligo_markers")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "00_check_Oligo_markers")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"))

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

sce <- sce[,sce$cell_type_broad %in% c("Oligo", "OPC")]

sce$cell_type_anno2 <- sce$cell_type_anno ## keep
sce$cell_type_anno <- as.character(sce$cell_type_anno)

## flatten OPC subtypes
sce$cell_type_anno[sce$cell_type_broad == "OPC"] <- "OPC"

sce$cell_type_anno <- factor(sce$cell_type_anno)
table(sce$cell_type_anno)

reducedDims <- readRDS(here("processed-data", "04_snRNA-seq", "38_sn_subcluster_reducedDims_OligoOPC", "Oligo_OPC_reducedDims.Rds"))
names(reducedDims)

reducedDim(sce, "UMAP") <- reducedDims$UMAP 
reducedDim(sce, "TSNE") <- reducedDims$TSNE 

rownames(sce) <- rowData(sce)$gene_name

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)
# Oligo_OPC_colors
Oligo_OPC_colors2 <- Oligo_OPC_colors
Oligo_OPC_colors <- c(Oligo_OPC_colors[grepl("Oligo", names(Oligo_OPC_colors))], c(OPC = "#D2B037"))


#### Data Driven Markers ####
list.files(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo"))
load(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "marker_stats_MeanRatio_Oligo.Rdata"), verbose = TRUE)
# marker_stats_MeanRatio

marker_stats_MeanRatio |>
    filter(MeanRatio > 1, MeanRatio.rank <= 5)

enrichment_top100 <- read_csv(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "subtype_enrichment_top100_OligoOPC.csv"))

enrichment_top100 |> 
    group_by(cluster = test) |>
    arrange(-stat) |>
    slice(1:5)

#### Violin plots ####

lit_marker_sets <- list(Mathys2024 = c("OPALIN", "RASGRF1"),
                        Siletti2023 = c("OPALIN", "RBFOX1", "CD44"))

plot_marker_express_List(sce,
                        gene_list = lit_marker_sets,
                        cellType_col = "cell_type_anno",
                        color_pal = Oligo_OPC_colors,
                        pdf_fn = here(plot_dir, "Violin_Oligo_lit_markers.pdf")
                        )

#### Plot markers in reduced Dims ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

walk(c("OPALIN", "LINGO2", "PLP1", "RASGRF1","RBFOX1", "CD44", "ARHGEF3"),
     ~my_plot_reduced_dim(sce,
                          prefix = "sn_subcluster",
                          dimred = "TSNE",
                          my_var = .x,
                          var_type = "express",
                          save_plot = TRUE,
                          facet = FALSE,
                          plot_dir_rd = plot_dir,
                          verbose = TRUE,
                          NA_gray = TRUE)
)



