## Louise Huuki-Myers, Match 2025
## Check Oligo marker genes (data drven or lit)

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("DeconvoBuddies")
library("scDotPlot")
library("readxl")

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

sce$cell_type_anno2 <- factor(sce$cell_type_anno2, levels = c("OPC.3", "OPC.4", "OPC.1", "OPC.2", "OPC.5", "Oligo.3", "Oligo.5", "Oligo.4", "Oligo.1", "Oligo.2"))

sce$cell_type_anno <- as.character(sce$cell_type_anno)

## flatten OPC subtypes
sce$cell_type_anno[sce$cell_type_broad == "OPC"] <- "OPC"

sce$cell_type_anno <- factor(sce$cell_type_anno, levels = c("OPC", "Oligo.3", "Oligo.5", "Oligo.4", "Oligo.1", "Oligo.2"))
table(sce$cell_type_anno)

table(sce$cell_type_anno, sce$cell_type_anno2)

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
# Oligo v Oligo
list.files(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo"))

## MeanRatio
load(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "marker_stats_MeanRatio_Oligo.Rdata"), verbose = TRUE)

marker_stats_MeanRatio |> 
    filter(gene %in% c("RASGRF1", "RBFOX1", "CD44"))

MeanRatio_Oligo <- marker_stats_MeanRatio |>
    filter(MeanRatio > 1, MeanRatio.rank <= 10)  |>
    select(cellType = cellType.target, gene) |>
    mutate(method = "MeanRatio_Oligo")

## enrich
enrich_Oligo <- read_csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_enrichment_top100_Oligo.csv")) 

enrich_Oligo|>
    filter(gene %in% c("RASGRF1", "RBFOX1", "CD44"))

enrich_Oligo_filter <- enrich_Oligo |>
    filter(top <= 10) |>
    select(cellType = test, gene) |>
    mutate(method = "Enrich_Oligo")


## Oligo OPC
# list.files(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC"))
## MeanRatio
load(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "marker_stats_MeanRatio_OligoOPC.Rdata"), verbose = TRUE)

marker_stats_MeanRatio |> 
    filter(gene %in% c("RASGRF1", "RBFOX1", "CD44"))

MeanRatio_OligoOPC <- marker_stats_MeanRatio |>
    filter(MeanRatio > 1, MeanRatio.rank <= 10)  |>
    select(cellType = cellType.target, gene) |>
    mutate(method = "MeanRatio_OligoOPC")

## enrich
enrich_OligoOPC <- read_csv(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "subtype_enrichment_top100_OligoOPC.csv"))

enrich_OligoOPC|>
    filter(gene %in% c("RASGRF1", "RBFOX1", "CD44"))

enrich_OligoOPC_filter <- enrich_OligoOPC |>
    filter(top <= 10) |>
    select(cellType = test, gene) |>
    mutate(method = "Enrich_OligoOPC")


Oligo_markers_long <- MeanRatio_Oligo |>
    bind_rows(enrich_Oligo_filter) |>
    bind_rows(MeanRatio_OligoOPC) |>
    bind_rows(enrich_OligoOPC_filter) |>
    group_by(cellType, gene) |>
    summarise(n = n(),
              methods = paste(method, collapse = ","))

Oligo_markers_long_2plus <- Oligo_markers_long |>  
    filter(n > 1) 

Oligo_markers_long_2plus |> count(cellType)


Oligo_markers_long |> arrange(-n) |> filter(cellType == "Oligo.3") |> print(n = 30)

Oligo_markers_long |> count(cellType, methods) |> print(n )

Oligo_markers_long |> ungroup() |> group_by(gene) |> filter(n() > 1)

Oligo_marker_summary <- MeanRatio_Oligo |>
    group_by(cellType) |>
    summarise(MeanRatio_Oligo = paste(gene, collapse = ","))|> 
    left_join(enrich_Oligo_filter |>
                  group_by(cellType) |>
                  summarise(Enrichment_Oligo = paste(gene, collapse = ","))) |> 
    full_join(MeanRatio_OligoOPC |>
                  group_by(cellType) |>
                  summarise(MeanRatio_OligoOPC = paste(gene, collapse = ","))) |>
    left_join(enrich_OligoOPC_filter |>
                  group_by(cellType) |>
                  summarise(Enrichment_OligoOPC = paste(gene, collapse = ","))) |>
    left_join(Oligo_markers_long_2plus |>
                  group_by(cellType) |>
                  summarise(multiHit = paste(gene, collapse = ",")))

write_csv(Oligo_marker_summary, file = here(data_dir, "ERC_Oligo_marker_summary.csv"))

#### Violin plots ####
## TODO fininsh annotations
erc_oligo_key_genes <- list(OPC = c("PDGFRA", "MEG3","OLIG2"), #OPCs
                         "RBFOX1", "KCND2", "GPM6A", #OPC + Oligo.3
                         "CNTNAP2", "NTRK3","KCNJ3", #OPC + Oligo.3
                         "LINGO2", "MT-CO3", "ADGRV1",   # Oligo.3
                         "ARHGEF3", "ADGRF5",  "CLDN5", #Oligo.5
                         "OPALIN", "OMG", "SEMA6D", #Oligo.1
                         "LAMA2", "ERBB4",  #Oligo.1 + 4
                         "RASGRF1","RASGRF2", "LRRC63", "ANKRD18A" #Oligo.2
)


violin_key_genes <- plot_gene_express(sce, 
                  category = "cell_type_anno",
                  genes = erc_oligo_key_genes,
                  color_pal = Oligo_OPC_colors,
                  free_y = TRUE,
                  ncol = 3
                  )

ggsave(violin_key_genes, filename = here(plot_dir, "Violin_Oligo_key_markers.png"))

violin_OPLAIN_RASGRF1_genes <- plot_gene_express(sce, 
                                      category = "cell_type_anno",
                                      genes = c("OPALIN", "RBFOX1", "RASGRF1"),
                                      color_pal = Oligo_OPC_colors,
                                      free_y = TRUE,
                                      ncol = 3
)

ggsave(violin_OPLAIN_RASGRF1_genes, filename = here(plot_dir, "Violin_Oligo_OPLAIN_RASGRF1.png"), width = 8, height = 4)


rowData(sce)$cop_marker <- NULL
rowData(sce)$cop_marker <- names(fang_cop_markers2)[match(rownames(sce), fang_cop_markers2)] 
table(rowData(sce)$cop_marker)

pdf(here(plot_dir, "OligoOPC_dotplot_key_markers.pdf"))
sce |>
    scDotPlot(features = erc_oligo_key_genes,
              group = "cell_type_anno2",
              groupAnno = "cell_type_anno2",
              # featureAnno = "cop_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno2" = Oligo_OPC_colors2),
              clusterRows = TRUE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = erc_oligo_key_genes,
              group = "cell_type_anno2",
              groupAnno = "cell_type_anno2",
              # featureAnno = "cop_marker",
              scale = FALSE,
              annoColors = list("cell_type_anno2" = Oligo_OPC_colors2),
              clusterRows = TRUE,
              groupLegends = FALSE)

dev.off()

pdf(here(plot_dir, "Oligo_dotplot_key_markers.pdf"))
sce[,sce$cell_type_broad == "Oligo"] |>
    scDotPlot(features = erc_oligo_key_genes,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "cop_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = Oligo_OPC_colors),
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()


## Marques
Marques_markers <- consensus_marker_sets <- list(
    OPC = c("PTPRZ1", "PDGFRA", "SERPINE2", "CSPG4", "CSPG5", "VCAN"),
    COP = c("CD9", "NEU4", "BMP4","GPR17","VCAN"),
    NFOL = c("ARPC1B", "TMEM2", "CHN2", "MPZL1", "FRMD4A", "MOBP", "DDR1", "TSPAN2"),
    MFOL = c("CTPS1", "TMEM141", "OPALIN", "MAL", "PTGDS", "EVI2A", "EVI2B"),
    MOL = c("APOD", "SEPP1", "S100B"),
    MOL1 = c("FOSB", "DUSP1", "DNAJB1"),
    MOL2 = c("ANXA5", "KLK6", "MGST3"),
    MOL3 = c("CAR2", "CNTN2", "GAD2"),
    MOL4 = c("SERPINB1", "NEAT1"),
    MOL5 = c("CYP51A1", "DHCR24", "PDLIM2"),
    MOL6 = c("IL33", "APOE", "PTGDS")
)


lit_marker_sets <- list(Marques2016 = c("PDGFRA", "CSPG4", "PTPRZ1", "PCDH15", "VCAN", "SOX6", "GPR17", "NEU4", "BMP4", "NKX2-2"),
                        Jakel2019 = c("SOX6", "BCAN", "ITPR2", "APOE", "CD74", "KLK6", "OPALIN", "PLP1"), 
                        Mathys2024 = c("OPALIN", "RASGRF1"),
                        Siletti2023 = c("OPALIN", "RBFOX1", "CD44", "GPC5"))

plot_marker_express_List(sce,
                        gene_list = lit_marker_sets,
                        cellType_col = "cell_type_anno",
                        color_pal = Oligo_OPC_colors,
                        pdf_fn = here(plot_dir, "Violin_Oligo_lit_markers.pdf")
                        )

#### Oligo + OPC subtypes #####

violin_Siletti2023 <- plot_gene_express(sce, 
                                      category = "cell_type_anno2",
                                      genes = c("HES5", "IRX5", "GPC5", "OPALIN", "RBFOX1", "CD44", "GPC5"),
                                      color_pal = Oligo_OPC_colors2,
                                      free_y = TRUE
)

ggsave(violin_Siletti2023, filename = here(plot_dir, "Violin_Oligo_Siletti2023_markers.png"))


pdf(here(plot_dir, "Oligo_dotplot_Siletti2023.pdf"))
sce |>
    scDotPlot(features = c("HES5", "IRX5", "GPC5", "OPALIN", "RBFOX1", "RASGRF1", "CD44", "LINGO2"),
              group = "cell_type_anno2",
              groupAnno = "cell_type_anno2",
              # featureAnno = "Marques_markers",
              scale = TRUE,
              annoColors = list("cell_type_anno2" = Oligo_OPC_colors2),
              clusterRows = TRUE,
              groupLegends = FALSE)
dev.off()


#### Plot markers in reduced Dims ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

# walk(c("OPALIN", "LINGO2", "PLP1", "RASGRF1","RBFOX1", "CD44", "ARHGEF3"),
walk(c("GPR17", "BCAS1"),
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


#### COP ####
# Fang et al., 2023 10.1002/glia.24426

fang_cop_markers <- list(OPC = c("NG2", "PDGFRA"),
                         COP = c("GPR17", "BCAS1", "FYN", "NKX2.2", "O4", 
                                 "EPHB1", "SH3RF3" #10.1016/j.neuron.2022.09.010
                                 ),
                         NFO = c("MYRF", "PLP", "CNP", "MBP", "MAG", "CC1"),
                         MOL = c("ASPA", "TMEM10"))

fang_cop_markers2 <- map(fang_cop_markers, ~.x[.x %in% rownames(sce)])

plot_marker_express_List(sce,
                         gene_list = fang_cop_markers2,
                         cellType_col = "cell_type_anno2",
                         color_pal = Oligo_OPC_colors2,
                         pdf_fn = here(plot_dir, "Violin_OligoOPC_fang_cop_markers.pdf")
)

## oligo marker dot plot
fang_cop_markers2 <- AnnotationDbi::unlist2(fang_cop_markers2)

rowData(sce)$cop_marker <- NULL
rowData(sce)$cop_marker <- names(fang_cop_markers2)[match(rownames(sce), fang_cop_markers2)] 
table(rowData(sce)$cop_marker)

pdf(here(plot_dir, "Oligo_dotplot_COP_Fang2023.pdf"))
sce |>
    scDotPlot(features = fang_cop_markers2,
              group = "cell_type_anno2",
              groupAnno = "cell_type_anno2",
              featureAnno = "cop_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno2" = Oligo_OPC_colors2),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = fang_cop_markers2,
              group = "cell_type_anno2",
              groupAnno = "cell_type_anno2",
              featureAnno = "cop_marker",
              scale = FALSE,
              annoColors = list("cell_type_anno2" = Oligo_OPC_colors2),
              clusterRows = FALSE,
              groupLegends = FALSE)

dev.off()

#### Siletti2023 z-stat ####

Siletti2023_markers <- read_xlsx(here("external-data", "Siletti2023", "science.add7046_table_s3.xlsx")) |>
    filter(grepl("Oligo|oligo", Supercluster)) |>
    select(Supercluster, `Cluster name`, `Number of cells`, `Top Enriched Genes`) |>
    mutate(Genes = str_split(`Top Enriched Genes`, pattern = ", "),
           n_genes = length(Genes))

# message(Sys.time(), "- Pseudobulk OligoOPC data")
# sce_pb <- scuttle::aggregateAcrossCells(
#     sce,
#     DataFrame(
#         cll_type_broad = sce[["cell_type_broad"]],
#         cell_type_anno = sce[["cell_type_anno"]],
#         cell_type_anno2 = sce[["cell_type_anno2"]]
#     )
# )
# 
# colnames(sce_pb) <- sce_pb$cell_type_anno2
# logcounts(sce_pb)  <-
#     edgeR::cpm(
#         edgeR::calcNormFactors(sce_pb),
#         log = TRUE,
#         prior.count = 1
#     )

# saveRDS(sce_pb, file = here(data_dir, "sce_pb_OligoOPC.rds"))
# message(Sys.time(), "- Done")

sce_pb <- readRDS(file = here(data_dir, "sce_pb_OligoOPC.rds"))


Siletti2023_markers_list <- Siletti2023_markers$Genes
names(Siletti2023_markers_list) <- Siletti2023_markers$`Cluster name`

Siletti2023_markers_list <- map(Siletti2023_markers_list, ~.x[.x%in% rownames(sce)])
map_int(Siletti2023_markers_list, length)

logcounts(sce_pb)[Siletti2023_markers_list$COP_75,]

spatialLIBD:::multi_gene_z_score(t(logcounts(sce_pb)[Siletti2023_markers_list$COP_37,]))

Siletti_z_score <- map_dfr(Siletti2023_markers_list, ~spatialLIBD:::multi_gene_z_score(t(logcounts(sce_pb)[.x,])))

Siletti_z_score <- as.matrix(Siletti_z_score)
rownames(Siletti_z_score) <- names(Siletti2023_markers_list)

pdf(here(plot_dir, "OligoOPC_Siletti2023_marker_z-score.pdf"), height = 4, width = 8)
ComplexHeatmap::Heatmap(Siletti_z_score, name = "marker z-score")
dev.off()

#### gene expression vs. pseudotime ####

oligo_lineage_markers <- tibble(
    gene = c("PDGFRA", "CSPG4", "OLIG2", "SOX10", "CNP", "MBP", "PLP1", "PTGSD", "RBFOX1", "RASGRF1", "OMG", "LINGO2", "GPM6A", "OPALIN", "MOG"),
    target = c("OPC", "OPC", "OPC", "Oligo.3", "Oligo", "Oligo", "Oligo", "Oligo", "OPC", "Oligo", "Oligo", "Oligo.3", "OPC", "Oligo", "Oligo")) |> 
  filter(gene %in% rownames(sce_pb))

lineage_z_score <- scale(t(logcounts(sce_pb)[oligo_lineage_markers$gene,]))


lineage_z_score_long <- lineage_z_score |>
  as.data.frame() |>
  rownames_to_column("cell_type") |>
  pivot_longer(!cell_type, names_to = "gene", values_to = "z_score") |>
  mutate(cell_type = factor(cell_type, levels = levels(sce_pb$cell_type_anno2))) |>
  left_join(oligo_lineage_markers)

lineage_z_score_plot <- lineage_z_score_long |>
  ggplot(aes(x = cell_type, y = z_score, color = gene)) +
  geom_point() +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  geom_line(aes(group = gene)) +
  facet_wrap(~target, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(lineage_z_score_plot, filename = here(plot_dir, "lineage_z_score_line.png"))

readRDS(here("processed-data", "04_snRNA-seq", "38_sn_subcluster_reducedDims_OligoOPC", "Oligo_OPC_line.data.Rds"))

#### Complie other dataset correlatons ####

dataset_cor_fn  <- map(c(Sadick2022 =  "sadick", Grubman2019 = "grubman"), ~here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", sprintf("erc_v_%s_oligo_cor.csv", .x)))

## jakel
dataset_cor_fn$Jakel2019 <- here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "erc_v_jakel_OligoOPC_cor.csv")

dataset_cor_fn$Green2024 <- here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC2","erc_v_Green2024_OligoOPC2_cor.csv") 

dataset_cor <- map2(dataset_cor_fn, names(dataset_cor_fn), ~read_csv(.x) |> mutate(dataset = .y))
dataset_cor$Sadick2022  <- dataset_cor$Sadick2022 |> dplyr::filter(data == "Sadick") |> select(-data) |> mutate(Sadick_cluster = gsub("Sadick_", "", Sadick_cluster))
dataset_cor$Jakel2019 <- dataset_cor$Jakel2019 |> mutate(Jakel_Oligo = gsub("Jakel_", "", Jakel_Oligo))

dataset_cor_all <- map2_dfr(dataset_cor, names(dataset_cor), ~.x |> 
  rename_with(~"other_oligo", .cols = 2))

write_csv(dataset_cor_all, file = here(data_dir, "ERC_Oligo_v_external_data_cor.csv"))


max_cor_oligo3 <-  dataset_cor_all |> 
  filter(test == "Oligo.3") |> 
  mutate(dataset = .y) |>
  slice_max(cor)) |>
  mutate(other_oligo = fct_reorder(paste0(dataset, ": ", other_oligo), cor)) |>
  arrange(other_oligo)

max_cor_oligo.3_barplot <- max_cor_oligo3 |>
  ggplot(aes(x = cor, y = other_oligo, fill = cor)) +
  geom_col() +
  theme_bw() +
  theme(legend.position = "None") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0,1)) +
  labs(x = "correlation\nlogFC marker stats", y = NULL) 

ggsave(max_cor_oligo.3_barplot, filename = here(plot_dir, "max_cor_oligo.3_barplot.png"), height = 3 , width = 3)
  
