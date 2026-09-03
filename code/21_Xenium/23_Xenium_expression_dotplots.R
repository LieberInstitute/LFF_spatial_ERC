## Louise Huuki-Myers, July 2025
## Make dotplots of select gene sets

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("scDotPlot")
library("tidyverse")

plot_dir <- here("plots", "05_spe_correct_cluster", "30_SpD_dotplot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "30_SpD_dotplot")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

rownames(spe) <- rowData(spe)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
# load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
SpD_colors <- metadata(spe)$SpD_colors

spd_levels <- levels(spe$vSpD)

#### Sex gene dot plot ####
sex_check <- list(male = c("SRY", "RPS4Y1", "RPS4Y2", "DDX3Y", "KDM5D", "UTY", "ZFY", "EIF1AY", "USP9Y", "TSPY1"),
                  female = c("XIST", "TSIX", "KDM6A", "EIF2S3X", "RPS4X"))

sex_check <- map(sex_check, ~.x[.x %in% rownames(spe)])

sex_check_tb <- tibble(gene = unlist(sex_check), 
                       sex_check = unlist(map2(names(sex_check), sex_check, ~rep(.x, length(.y)))))

# Autosomal, Hormone-Responsive
# "FOXP3", "ESR1", "AR", "TSHR", "PRLR"

rowData(spe)$sex_check <- sex_check_tb$sex_check[match(rownames(spe), sex_check_tb$gene)] 
table(rowData(spe)$sex_check)


dotplot_sex_genes <- spe |>
    scDotPlot(features = unlist(sex_check),
              group = "sample_id",
              groupAnno = "Sex",
              featureAnno = "sex_check",
              scale = TRUE,
              annoColors = list("Sex" = sex_colors),
              clusterRows = FALSE,
              groupLegends = FALSE
              )

ggsave(dotplot_sex_genes, filename = here(plot_dir, "Visium_sample_dotplot_Sex_genes.pdf"))

#### Lit Dotplots ####

lit_markers <- read_csv(here("processed-data","05_spe_correct_cluster", "00_lit_marker_genes_layer", "lit_layer_marker_summary.csv")) |>
    filter(gene_name %in% rowData(spe)$gene_name) |>
    arrange(Layer)

rowData(spe)$LitMarker <- NULL
rowData(spe)$LitMarker <- lit_markers$Layer[match(rownames(spe), lit_markers$gene_name)] 
table(rowData(spe)$LitMarker)

layer_colors <- c(Vasc = '#FF99C0',
                  Layer1 = "#F0027F",
                  Layer2 = "#377EB8",
                  Layer3 = "#4DAF4A",
                  Layer4 = "#984EA3",
                  Layer5 = "#FFD700",
                  Layer5a = "#E9FF70",
                  Layer5b = "#FF9770",
                  Layer6 = "#FF7F00",
                  WM = "grey90",
                  GM = "grey25",
                  `Para-subiculum` = "brown",
                  Interneuron = "red")

dotplot_lit_markers <- spe |>
    scDotPlot(features = lit_markers$gene_name,
              group = "vSpD",
              groupAnno = "vSpD",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("vSpD" = SpD_colors,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_lit_markers, filename = here(plot_dir, "Visium_SpD_dotplot_LitMarkers.pdf"))
ggsave(dotplot_lit_markers, filename = here(plot_dir, "Visium_SpD_dotplot_LitMarkers.png"))

key_ramsden_markers <- c("NXPH4", "FEZF2", "ETV1", "DCC", "RELN", "WFS1")

dotplot_lit_markers_ramsden <- spe |>
    scDotPlot(features = key_ramsden_markers,
              group = "vSpD",
              groupAnno = "vSpD",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("vSpD" = SpD_colors,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_lit_markers_ramsden, filename = here(plot_dir, "Visium_SpD_dotplot_LitMarkers_Ramsden.pdf"), width = 4, height = 6)

dotplot_lit_markers_ramsden_all <- spe |>
    scDotPlot(features = lit_markers |> filter(grepl("Ramsden", studies)) |> pull(gene_name),
              group = "vSpD",
              groupAnno = "vSpD",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("vSpD" = SpD_colors,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_lit_markers_ramsden_all, filename = here(plot_dir, "Visium_SpD_dotplot_LitMarkers_Ramsden_all.pdf"), width = 7, height = 7)



## Ramsden

Ramsden_anno = c("L2", "L3", "LVa", "L5b", "L6", "Island")
Gene = c("DCC", "KITL", "ETV1", "FEZF2", "NXPH4", "WFS1")

Gene %in% rownames(spe)


#### MeanRatio dot plots ####
load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_vSpD.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(MeanRatio.rank <= 6, MeanRatio >1) |>
    select(gene, MeanRatio.rank, vSpD = cellType.target) |>
    mutate(vSpD = factor(vSpD, levels = spd_levels)) |>
    arrange(vSpD)

marker_stats_top |> count(vSpD)

rowData(spe)$vSpD_Marker <- NULL
rowData(spe)$vSpD_Marker <- marker_stats_top$vSpD[match(rownames(spe), marker_stats_top$gene)] 
table(rowData(spe)$vSpD_Marker)

dotplot_mean_ratio <- spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "vSpD",
              groupAnno = "vSpD",
              featureAnno = "vSpD_Marker",
              scale = TRUE,
              annoColors = list("vSpD" = SpD_colors,
                                "vSpD_Marker" = SpD_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_mean_ratio, filename = here(plot_dir, "Visium_SpD_dotplot_MeanRatio.pdf"), width = 7, height = 8)
ggsave(dotplot_mean_ratio, filename = here(plot_dir, "Visium_SpD_dotplot_MeanRatio.png"), width = 7, height = 8)


#### Enrichment Data ####
erc_modeling_top <- read_csv(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "All_modeling_SpD_top100.csv"))

erc_modeling_top |> filter(test == "vL5a-vL5b", fdr < 0.01) |> arrange(-logFC) |> print(n = 30)
erc_modeling_top |> filter(test == "vL5b-vL5a", fdr < 0.01) |> arrange(-logFC)

enrichment_top <- erc_modeling_top |>
    filter(model_type == "enrichment", top <=5) |>
    select(gene, top, vSpD = test) |>
    mutate(vSpD = factor(vSpD, levels = spd_levels)) |>
    arrange(vSpD)

enrichment_top |> count(vSpD)

rowData(spe)$vSpD_enrich <- NULL
rowData(spe)$vSpD_enrich <- enrichment_top$vSpD[match(rownames(spe), enrichment_top$gene)] 
table(rowData(spe)$vSpD_enrich)

dotplot_enrichment <- spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "vSpD",
              groupAnno = "vSpD",
              featureAnno = "vSpD_enrich",
              scale = TRUE,
              annoColors = list("vSpD" = SpD_colors,
                                "vSpD_enrich" = SpD_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_enrichment, filename = here(plot_dir, "Visium_SpD_dotplot_enrichment.pdf"), width = 7, height = 8)
ggsave(dotplot_enrichment, filename = here(plot_dir, "Visium_SpD_dotplot_enrichment.png"), width = 7, height = 8)

#### Key genes ####

erc_SpD_key <- read_csv(here(data_dir, "vSpD_key_genes_long.csv"))

rowData(spe)$vSpD_genes <- NULL
rowData(spe)$vSpD_genes <- erc_SpD_key$domain[match(rownames(spe), erc_SpD_key$gene)] 
table(rowData(spe)$vSpD_genes)

dotplot_key_genes <- spe |>
    scDotPlot(features = erc_SpD_key$gene,
              group = "vSpD",
              groupAnno = "vSpD",
              featureAnno = "vSpD_genes",
              scale = TRUE,
              annoColors = list("vSpD" = SpD_colors,
                                "vSpD_genes" = c(SpD_colors, lit = "black")),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_key_genes, filename = here(plot_dir, "Visium_SpD_dotplot_key_genes.pdf"), width = 5, height = 7)
ggsave(dotplot_key_genes, filename = here(plot_dir, "Visium_SpD_dotplot_key_genes.png"), width = 5, height = 7)


# slurmjobs::job_single('30_SpD_dotplot', create_shell = TRUE, memory = '5G', command = "Rscript 30_SpD_dotplot.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
