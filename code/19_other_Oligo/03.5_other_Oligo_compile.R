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
    select(ERC_Oligo = `Oligo subtype`, Oligo_anno = Annotation_short)

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)
Oligo_OPC_colors <- c(Oligo_OPC_colors[grepl("Oligo", names(Oligo_OPC_colors))], c(OPC = "#D2B037"))

Oligo_anno_colors <- c(Oligo.M  = "#00C4B8", 
                       Oligo.L  = "#0072E5",  
                       Oligo.NF = "#7B61FF", 
                       Oligo.P  = "#FF3D9A",  
                       OPC      = "#C49A00",
                       Oligo.other = "grey")


cor_layer_anno2_select_long <- cor_layer_anno_select_long |> 
    mutate(ERC_Oligo = gsub("/.*", "", layer_label)) |> 
    left_join(oligo_notes) |>
    replace_na(list(Oligo_anno = "Oligo.other")) |>
    mutate(dataset = gsub("_O.*","", cluster), .before = 1)


cor_layer_anno2_select_long |> dplyr::count(dataset, Oligo_anno)

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
    select(cluster, Oligo_anno)  |>
    column_to_rownames("cluster")

other_oligo_anno <- other_oligo_anno[rownames(cor_wide_matrix),]
other_oligo_ra <- rowAnnotation(" " = other_oligo_anno, col = list(` ` = Oligo_anno_colors), show_legend = FALSE)

identical(oligo_notes$ERC_Oligo, colnames(cor_wide_matrix))
      
erc_oligo_ca <- HeatmapAnnotation(Oligo_anno = oligo_notes$Oligo_anno, col = list(Oligo_anno = Oligo_anno_colors))

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

#### Top match heatmap ####
 
best_match <- cor_layer_anno2_select_long |> filter(ERC_Oligo == "Oligo.3")


cor_wide_matrix_best <- cor_wide_matrix[best_match$cluster,]

anno_matrix_best <- anno_matrix[best_match$cluster,]

pdf(here(plot_dir, "LIBD_Other_dataset_cor_bestOligo3.pdf"), height = 2.5, width = 4)
print(Heatmap(cor_wide_matrix_best, 
             name = "t-stat cor",
             col = circlize::colorRamp2(
                 breaks = seq(-1, 1, length.out = 7),
                 colors = RColorBrewer::brewer.pal(7, "PRGn")
             ),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(anno_matrix_best[i, j], x, y, gp = gpar(fontsize = 10))
             })
      )
dev.off()


#### Explore proportions ####

cluster_tab_fn <- list(spatialdACC = "walktrap_snn_k20_subclusters_spatialdACC_wOPC.Rdata",
                       spatialDLPFC = "walktrap_snn_k20_subclusters_spatialDLPFC_wOPC.Rdata",
                       spatialHPC = "walktrap_snn_k20_subclusters_spatialHPC.Rdata")

cluster_tab <- map_dfr(cluster_tab_fn, ~get(load(here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster", .x)))) |> 
    mutate(cluster = gsub("OPC.[0-9]","OPC", cluster_anno))

other_Oligo_cluster_count <- cor_layer_anno2_select_long |> 
    mutate(cell_type_broad = gsub(".*_|\\..*", "", cluster)) |> 
    select(dataset, cell_type_broad, cluster, Oligo_anno) |>
    full_join(cluster_tab |> 
                  dplyr::count(cluster))


cor_layer_anno2_select_long |> select(cluster)

erc_oligo_cluster_count <- read.csv(here("processed-data", "04_snRNA-seq", "33_sn_subcluster_summary", "ERC_sn_subcluster_summary_cell_type.csv"), row.names = 1) |>
    filter(cell_type_broad %in% c("Oligo", "OPC")) |>
    select(cell_type_broad, oligo_anno = cell_type_anno, n = n_nuclei) |>
    mutate(oligo_anno = gsub("OPC.*", "OPC", oligo_anno)) |>
    group_by(cell_type_broad, oligo_anno) |>
    summarise(n = sum(n)) |>
    group_by(cell_type_broad) |>
    mutate(oligo_anno = gsub("OPC.*", "OPC", oligo_anno), ## flatten OPCs
           prop = n/sum(n), 
           dataset= "erc", 
           ERC_Oligo = oligo_anno) |>
    left_join(oligo_notes)

other_erc_Oligo_cluster_count <- other_Oligo_cluster_count |>
    bind_rows(erc_oligo_cluster_count |> 
                  mutate(cluster = paste0("erc_", oligo_anno)) |>
                  select(dataset, cluster, Oligo_anno, n)) |>
    mutate(Oligo_anno = factor(Oligo_anno, levels = c("OPC", "Oligo.NF", "Oligo.P","Oligo.M", "Oligo.L", "Oligo.other")),
           level = as.integer(Oligo_anno),
           cluster = fct_reorder(cluster, level),
           dataset = factor(dataset, levels = c("erc", "dacc", "dlpfc", "hpc")))


Oligo_cluster_region_prop <- other_erc_Oligo_cluster_count |>
    group_by(dataset, cell_type_broad, Oligo_anno) |>
    summarize(n = sum(n))  |>
    group_by(dataset, cell_type_broad) |>
    mutate(prop = n/sum(n))
    

write_csv(Oligo_cluster_region_prop, here(data_dir, "Oligo_cluster_region_prop.csv"))

#### plot proportions ####

oligo_cluster_region_barplot <- other_erc_Oligo_cluster_count |>
    ggplot(aes(x = cluster, y = n, fill = Oligo_anno)) +
    geom_col() +
    facet_wrap(~dataset, scales = "free_x", nrow = 1, space = "free_x") +
    scale_fill_manual(values = Oligo_anno_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "n nuclei")

ggsave(oligo_cluster_region_barplot, filename = here(plot_dir, "oligo_cluster_region_n_barplot.png"), height = 3, width = 7)


oligo_cluster_region_prop_barplot <- Oligo_cluster_region_prop |>
    filter(cell_type_broad == "Oligo") |>
    ggplot(aes(x = dataset, y = prop, fill = Oligo_anno)) +
    geom_col() + 
    geom_text(aes(label = round(prop, 2)), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(Oligo_anno_colors)) +
    theme_bw() +
    labs(x = "Brain Region", y = "Proportion Oligo subtype") +
    coord_flip() +
    scale_x_discrete(limits = rev)

ggsave(oligo_cluster_region_prop_barplot, filename = here(plot_dir, "oligo_cluster_region_prop_barplot.png"), height = 2.5, width = 6)



    

         