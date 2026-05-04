## Louise Huuki-Myers, Dec 2025
## Plot gene expression of other Oligo clusters

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
# library("scDotPlot")

library("ComplexHeatmap")

data_dir <- here("processed-data", "19_other_Oligo", "03.5_other_Oligo_compile")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "03.5_other_Oligo_compile")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


cor_layer_fn <- list.files( here("processed-data", "19_other_Oligo", "02_other_Oligo_model"), pattern = "cor_layer_anno")

cor_layer_data <- map(cor_layer_fn, ~readRDS(here("processed-data", "19_other_Oligo", "02_other_Oligo_model", .x)))

# pdf(here(plot_dir, "layer_stat_cor_plot_all.pdf"))
# walk2(cor_layer_data, cor_layer_fn, ~print(layer_stat_cor_plot(cor_stats_layer = .x$cor_layer, annotation = .x$anno, column_title = basename(.y))))
# dev.off()


cor_layer_select_fn <-  list(
    HPC = "other_Oligo_cor_layer_anno_spatialHPC_k20.Rds",
    dACC = "other_OligoOPC_cor_layer_anno_spatialdACC_wOPC_k20.Rds",
    DLPFC = "other_OligoOPC_cor_layer_anno_spatialDLPFC_wOPC_k20.Rds")

cor_layer_data_select <- map(cor_layer_select_fn, ~readRDS(here("processed-data", "19_other_Oligo", "02_other_Oligo_model", .x)))


cor_layer_select_long <- map2_dfr(cor_layer_data_select, names(cor_layer_data_select), 
                                  ~.x$'cor_layer' |> 
                                      reshape2::melt() |> 
                                      as_tibble() |>
                                      dplyr::rename(Other_Oligo = Var1, ERC_Oligo = Var2, cor = value) |>
                                      mutate(dataset = .y,
                                             ERC_Oligo = factor(ERC_Oligo, 
                                                    levels = c("OPC", "Oligo.3", "Oligo.5", "Oligo.4", "Oligo.1", "Oligo.2")
                                             )))

cor_layer_anno_select_long <- map2_dfr(cor_layer_data_select, names(cor_layer_data_select), 
                                  ~.x$'anno')

#### other oligo annotation ####

oligo_notes <- readxl::read_xlsx(here("processed-data", "19_other_Oligo", "00_check_Oligo_markers", "ERC_Oligo_notes.xlsx")) |>
    select(ERC_Oligo = `Oligo subtype`, Annotation_short)

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)
Oligo_OPC_colors <- c(Oligo_OPC_colors[grepl("Oligo", names(Oligo_OPC_colors))], c(OPC = "#D2B037"))

Oligo_anno_colors <- c(Oligo.M = "#00BFC4",
                       # Oligo.M = "#E76BF3",
                       Oligo.L = "#00B0F6",
                       Oligo.NF = "#9590FF",
                       Oligo.P = "#FF62BC",
                       OPC = "#D2B037")

cor_layer_anno2_select_long <- cor_layer_anno_select_long |> 
    mutate(ERC_Oligo = gsub("/.*", "", layer_label)) |> 
    left_join(oligo_notes)

write_csv(cor_layer_anno2_select_long, file = here(data_dir, "LIBD_other_Oligo_annotations.csv"))



#### combined t-stat heatmp ####

cor_wide <- cor_layer_select_long |> 
    arrange(ERC_Oligo) |>
    select(dataset, ERC_Oligo, Other_Oligo, cor) |>
    pivot_wider(names_from = "ERC_Oligo", values_from = "cor")


cor_wide_matrix <- cor_wide |>
    select(-dataset) |>
    column_to_rownames("Other_Oligo") |>
    as.matrix()

## create X and * matrix 
anno_matrix <- spatialLIBD:::create_annotation_matrix(cor_layer_anno_select_long, cor_wide_matrix)

other_oligo_anno <- cor_layer_anno2_select_long |> 
    select(cluster, Oligo_anno = Annotation_short)  |>
    column_to_rownames("cluster")

other_oligo_anno <- other_oligo_anno[rownames(cor_wide_matrix),]
other_oligo_ra <- rowAnnotation(" " = other_oligo_anno, col = list(` ` = Oligo_anno_colors), show_legend = FALSE)

identical(oligo_notes$ERC_Oligo, colnames(cor_wide_matrix))
      
erc_oligo_ca <- HeatmapAnnotation(Oligo_anno = oligo_notes$Annotation_short, col = list(Oligo_anno = Oligo_anno_colors))


pdf(here(plot_dir, "LIBD_Other_dataset_cor.pdf"), height = 8, width = 5)
draw(Heatmap(cor_wide_matrix, 
              name = "t-stat cor",
              col = circlize::colorRamp2(
                  breaks = seq(-1, 1, length.out = 7),
                  colors = RColorBrewer::brewer.pal(7, "PRGn")
              ),
              right_annotation = other_oligo_ra,
              bottom_annotation = erc_oligo_ca,
              row_split = cor_wide$dataset,
              cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(anno_matrix[i, j], x, y, gp = gpar(fontsize = 10))
              }
), merge_legend = TRUE,)
dev.off()


         