## Louise Huuki-Myers, March 2025
## Compile Oligo marker stat correlation from external data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ComplexHeatmap")

data_dir <- here("processed-data", "19_other_Oligo", "00_other_Oligo_cor")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "00_other_Oligo_cor")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Complie other dataset correlatons ####

dataset_cor_fn  <- map(c(Sadick2022 =  "sadick", Grubman2019 = "grubman"), ~here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", sprintf("erc_v_%s_oligo_cor.csv", .x)))

## jakel
dataset_cor_fn$Jakel2019 <- here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "erc_v_jakel_OligoOPC_cor.csv")
dataset_cor_fn$Macnair2025 <- here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "erc_v_Macnair2025_OligoOPC_cor.csv")

dataset_cor_fn$Green2024 <- here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC","erc_v_Green2024_OligoOPC_cor.csv") 

dataset_cor <- map2(dataset_cor_fn, names(dataset_cor_fn), ~read_csv(.x) |> mutate(dataset = .y))
dataset_cor$Sadick2022  <- dataset_cor$Sadick2022 |> dplyr::filter(data == "Sadick") |> select(-data) |> mutate(Sadick_cluster = gsub("Sadick_", "", Sadick_cluster))
dataset_cor$Jakel2019 <- dataset_cor$Jakel2019 |> mutate(Jakel_Oligo = gsub("Jakel_", "", Jakel_Oligo))

dataset_cor_all <- map2_dfr(dataset_cor, names(dataset_cor), ~.x |> 
                              rename_with(~"Other_Oligo", .cols = 2)) |>
  dplyr::rename(ERC_Oligo = test, n_gene = n, cor_logFC= cor)

dataset_cor_all |> filter(dataset == "Macnair2025") |> count(dataset, Other_Oligo)

## specificity score

# other_oligo_match <- dataset_cor_all |>
#   filter(grepl("Oligo", ERC_Oligo)) |>
#   group_by(dataset, Other_Oligo) |>
#   mutate(rank_cor = rank(cor_logFC)) 
# 
# other_oligo_best_match <- other_oligo_match |>
#   slice_max(rank_cor) |>
#   dplyr::rename(best_cor_logFC = cor_logFC) |>
#   left_join(other_oligo_match |>
#               filter(rank_cor == 4) |>
#               select(-rank_cor, -n_gene, ERC_Oligo_next = ERC_Oligo) |>
#               dplyr::rename(next_cor_logFC = cor_logFC)) |>
#   mutate(specificity_score = (best_cor_logFC - next_cor_logFC)/abs(next_cor_logFC))
# 
# other_oligo_best_match|>
#   arrange(-specificity_score) |>
#   filter(dataset == "Sadick2022")
# 
# write_csv(other_oligo_best_match, file = here(data_dir, "ERC_Oligo_v_external_data_specificity_score.csv"))


# dataset_cor_all2 <- dataset_cor_all |>
#   left_join(other_oligo_best_match |> select(dataset, ERC_Oligo, Other_Oligo, specificity_score))

write_csv(dataset_cor_all, file = here(data_dir, "ERC_Oligo_v_external_data_cor.csv"))


# dataset_subtypes <- dataset_cor_all2 |> select(dataset, Other_Oligo) |> unique()
# write_csv(dataset_subtypes, file = here(data_dir, "Oligo_external_data_subtypes.csv"))

dataset_subtypes <- read_csv(here(data_dir, "Oligo_external_data_subtypes.csv"))
dataset_subtypes |> count(dataset)
dataset_subtypes |> count(marker)

# Oligo.1   Oligo.2   Oligo.3   Oligo.4   Oligo.5 
# #00BFC4" "#00B0F6" "#9590FF" "#E76BF3" "#FF62BC" 

Oligo_marker_colors <- c(OPALIN = "#228277",
                         RBFOX1 = "#05678e",
                         `RASGRF1/RBFOX1` = "#00B0F6",
                         RASGRF1 = "#9dd2e7")


#### uniform heatmap ####
dataset_groups <- list(Oligo = c("Sadick2022", "Grubman2019"), OligoOPC = c("Jakel2019", "Green2024", "Macnair2025"))

map2(dataset_groups, names(dataset_groups), function(datasets, name){
  
  other_data <- dataset_cor_all |> 
    filter(dataset %in% datasets) |>
    mutate(ERC_Oligo = droplevels(factor(ERC_Oligo, 
                                         levels = c("OPC.3", "OPC.4", "OPC.1", "OPC.2", "OPC.5", "OPC", "Oligo.3", "Oligo.5", "Oligo.4", "Oligo.1", "Oligo.2")
                                         )
                                  )
           )
  
  cor_wide <- other_data |> 
    arrange(ERC_Oligo) |>
    select(dataset, ERC_Oligo, Other_Oligo, cor_logFC) |>
    pivot_wider(names_from = "ERC_Oligo", values_from = "cor_logFC") |>
      mutate(Other_Oligo = ifelse(Other_Oligo == "COP" & dataset == "Green2024", "C0P", Other_Oligo)) ## double COP fix...
  
  # cor_wide |> count(Other_Oligo) |> filter(n>1)
  # cor_wide |> filter(Other_Oligo == "COP")
  
  cor_wide_matrix <- cor_wide |>
    select(-dataset) |>
    column_to_rownames("Other_Oligo") |>
    as.matrix()
  
  other_oligo_anno <- dataset_subtypes |> 
      filter(dataset %in% datasets) |>
      mutate(Other_Oligo = ifelse(Other_Oligo == "COP" & dataset == "Green2024", "C0P", Other_Oligo)) |>
      select(Other_Oligo, marker)  |>
      column_to_rownames("Other_Oligo")
  
  other_oligo_anno <- other_oligo_anno[rownames(cor_wide_matrix),]
  other_oligo_ra <- rowAnnotation(express = other_oligo_anno, col = list(express = Oligo_marker_colors))
  
  pdf(here(plot_dir, sprintf("Other_dataset_cor_%s.pdf", name)), height = 8, width = 5)
  print(Heatmap(cor_wide_matrix, 
                name = "logFC cor",
                right_annotation = other_oligo_ra,
                row_split = cor_wide$dataset,
                cluster_columns = FALSE
  ))
  dev.off()
  
})


max_cor_oligo3 <-  dataset_cor_all |> 
  filter(ERC_Oligo == "Oligo.3") |> 
  group_by(dataset) |>
  slice_max(cor_logFC)

# max_cor_oligo3 <-  dataset_cor_all |> 
#   filter(ERC_Oligo == "Oligo.3") |> 
#   group_by(dataset) |>
#   slice_max(cor_logFC) |> 
#   mutate(selection = "Max Cor") |>
#   bind_rows(dataset_cor_all2 |> 
#               filter(ERC_Oligo == "Oligo.3") |> 
#               group_by(dataset) |>
#               slice_max(specificity_score)|> 
#               mutate(selection = "Max Specificity") ) |>
#   group_by(Other_Oligo) |>
#   mutate(selection = ifelse(n() > 1, "Max Both", selection)) |>
#   unique() |>
#   mutate(Other_Oligo2 = paste0(dataset, ":", Other_Oligo)) |>
#   arrange(dataset)


#   mutate(Other_Oligo = fct_reorder(paste0(dataset, ": ", other_oligo), cor)) |>
#   arrange(other_oligo)

max_cor_oligo.3_barplot <- max_cor_oligo3 |>
  ggplot(aes(x = cor, y = Other_Oligo, fill = cor)) +
  geom_col() +
  theme_bw() +
  theme(legend.position = "None") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1)) +
  labs(x = "correlation\nlogFC marker stats", y = NULL) 

ggsave(max_cor_oligo.3_barplot, filename = here(plot_dir, "max_cor_oligo.3_barplot.png"), height = 3 , width = 3)

selection_colors = c(`Max Cor` = "#777acd",
                     `Max Both` = "#c45ca2",
                     `Max Specificity` = "#cb5a4c")

max_cor_oligo.3_dotplot <- max_cor_oligo3 |>
  filter(!is.na(specificity_score)) |>
  ggplot(aes(x = cor_logFC, y = Other_Oligo, size = specificity_score, color = `selection`)) +
  geom_point() +
  scale_color_manual(values = selection_colors) +
  theme_bw() +
  facet_grid(dataset~., scales = "free_y") +
  labs(x = "correlation\nlogFC marker stats", y = NULL) 

ggsave(max_cor_oligo.3_dotplot, filename = here(plot_dir, "max_cor_oligo.3_dotplot.png"), height = 5 , width = 3)


## combined heatmap 
oligo3_select_cor_matrix <- dataset_cor_all |>
  filter(Other_Oligo %in% max_cor_oligo3$Other_Oligo,
         grepl("Oligo", ERC_Oligo)) |>
    mutate(Other_Oligo = paste0(dataset, ":", Other_Oligo)) |>
  select(ERC_Oligo, Other_Oligo, cor_logFC) |>
  pivot_wider(values_from = cor_logFC, names_from = Other_Oligo) |>
  column_to_rownames("ERC_Oligo") |>
  as.matrix() |>
  t()

oligo3_select_cor_matrix <- oligo3_select_cor_matrix[order(oligo3_select_cor_matrix[,'Oligo.3']),]

setequal(rownames(oligo3_select_cor_matrix),max_cor_oligo3$Other_Oligo)
setdiff(rownames(oligo3_select_cor_matrix),max_cor_oligo3$Other_Oligo)

# oligo3_select_cor_matrix <- oligo3_select_cor_matrix[max_cor_oligo3$Other_Oligo,]
# 
# other_oligo_anno <- rowAnnotation(df = max_cor_oligo3 |> ungroup() |> select(selection)
# )

pdf(here(plot_dir, "Other_Oligo_select_heatmap.pdf"), height = 2.5, width = 4)
Heatmap(oligo3_select_cor_matrix,
        name = "cor logFC",
        cluster_rows = FALSE)
dev.off()

