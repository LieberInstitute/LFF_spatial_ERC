## Louise Huuki-Myers, July 2025
## Make dotplots of select gene sets
## Adapted from 05_spe_correct_cluster/30_SpD_dotplot.R for Xenium xSpD data

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("qs2")
library("scDotPlot")
library("tidyverse")

plot_dir <- here("plots", "21_Xenium", "23_Xenium_expression_dotplots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

data_dir <- here("processed-data", "21_Xenium", "23_Xenium_expression_dotplots")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load data ####
message(Sys.time(), " - Load Xenium SPE")
## xSpD version of spe, NOT filtered to singlets (matches 14_xenium_pseudobulk_model_register.R)
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "spe_xenium_bansky.qs2"))
spe

rownames(spe) <- rowData(spe)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
SpX_colors <- metadata(spe)$SpX_colors

spx_levels <- levels(spe$xSpD)

# #### Sex gene dot plot ####
# sex_check <- list(male = c("SRY", "RPS4Y1", "RPS4Y2", "DDX3Y", "KDM5D", "UTY", "ZFY", "EIF1AY", "USP9Y", "TSPY1"),
#                   female = c("XIST", "TSIX", "KDM6A", "EIF2S3X", "RPS4X"))
# 
# sex_check <- map(sex_check, ~.x[.x %in% rownames(spe)])
# 
# sex_check_tb <- tibble(gene = unlist(sex_check), 
#                        sex_check = unlist(map2(names(sex_check), sex_check, ~rep(.x, length(.y)))))
# 
# # Autosomal, Hormone-Responsive
# # "FOXP3", "ESR1", "AR", "TSHR", "PRLR"
# 
# rowData(spe)$sex_check <- sex_check_tb$sex_check[match(rownames(spe), sex_check_tb$gene)] 
# table(rowData(spe)$sex_check)
# 
# 
# dotplot_sex_genes <- spe |>
#     scDotPlot(features = unlist(sex_check),
#               group = "sample_id",
#               groupAnno = "Sex",
#               featureAnno = "sex_check",
#               scale = TRUE,
#               annoColors = list("Sex" = sex_colors),
#               clusterRows = FALSE,
#               groupLegends = FALSE
#               )
# 
# ggsave(dotplot_sex_genes, filename = here(plot_dir, "Xenium_sample_dotplot_Sex_genes.pdf"))

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
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("xSpD" = SpX_colors,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_lit_markers, filename = here(plot_dir, "Xenium_xSpD_dotplot_LitMarkers.pdf"))
ggsave(dotplot_lit_markers, filename = here(plot_dir, "Xenium_xSpD_dotplot_LitMarkers.png"))

## limit to genes present on the Xenium panel (custom, targeted panel - unlike whole-transcriptome Visium)
key_ramsden_markers <- c("NXPH4", "FEZF2", "ETV1", "DCC", "RELN", "WFS1")
key_ramsden_markers <- key_ramsden_markers[key_ramsden_markers %in% rownames(spe)]

dotplot_lit_markers_ramsden <- spe |>
    scDotPlot(features = key_ramsden_markers,
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("xSpD" = SpX_colors,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_lit_markers_ramsden, filename = here(plot_dir, "Xenium_xSpD_dotplot_LitMarkers_Ramsden.pdf"), width = 4, height = 6)

dotplot_lit_markers_ramsden_all <- spe |>
    scDotPlot(features = lit_markers |> filter(grepl("Ramsden", studies)) |> pull(gene_name),
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "LitMarker",
              scale = TRUE,
              annoColors = list("xSpD" = SpX_colors,
                                LitMarker = layer_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_lit_markers_ramsden_all, filename = here(plot_dir, "Xenium_xSpD_dotplot_LitMarkers_Ramsden_all.pdf"), width = 7, height = 7)



## Ramsden
Ramsden_anno = c("L2", "L3", "LVa", "L5b", "L6", "Island")
Gene = c("DCC", "KITL", "ETV1", "FEZF2", "NXPH4", "WFS1")

Gene %in% rownames(spe)


#### MeanRatio dot plots ####
## No Xenium-specific MeanRatio has been computed for xSpD; reuse the
## Visium vSpD MeanRatio markers (validated cross-platform), restricted to
## genes present on the Xenium panel - same approach used in
## 13.5_xenium_bansky_annotate.R and 15_xenium_vis.R. Genes are grouped on the
## x-axis by the Xenium xSpD domains, while featureAnno records which vSpD
## (Visium) domain each gene is a top marker for.
load(here("processed-data", "05_spe_correct_cluster", "27_SpD_MeanRatio", "marker_stats_MeanRatio_vSpD.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

marker_stats_top <- marker_stats |>
    filter(gene %in% rownames(spe), MeanRatio > 1) |>
    group_by(cellType.target) |>
    arrange(MeanRatio.rank) |>
    slice(1:6) |>
    ungroup() |>
    select(gene, MeanRatio.rank, vSpD = cellType.target) |>
    arrange(vSpD)

marker_stats_top |> count(vSpD)

rowData(spe)$vSpD_Marker <- NULL
rowData(spe)$vSpD_Marker <- marker_stats_top$vSpD[match(rownames(spe), marker_stats_top$gene)] 
table(rowData(spe)$vSpD_Marker)

dotplot_mean_ratio <- spe |>
    scDotPlot(features = marker_stats_top$gene,
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "vSpD_Marker",
              scale = TRUE,
              annoColors = list("xSpD" = SpX_colors,
                                "vSpD_Marker" = SpD_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_mean_ratio, filename = here(plot_dir, "Xenium_xSpD_dotplot_MeanRatio.pdf"), width = 7, height = 8)
ggsave(dotplot_mean_ratio, filename = here(plot_dir, "Xenium_xSpD_dotplot_MeanRatio.png"), width = 7, height = 8)


#### Enrichment Data ####
## Use enrichment modeling run natively on Xenium xSpD clusters
## (14_xenium_pseudobulk_model_register.R --var xSpD), so domain names/colors
## line up directly with spe$xSpD / SpX_colors.
erc_modeling_top <- read_csv(here("processed-data", "21_Xenium", "14_xenium_pseudobulk_model_register", "xenium_enrichment_modeling_xSpD_top100.csv"))

enrichment_top <- erc_modeling_top |>
    filter(model_type == "enrichment", top <=5) |>
    select(gene, top, xSpD = test) |>
    mutate(xSpD = factor(xSpD, levels = spx_levels)) |>
    arrange(xSpD)

enrichment_top |> count(xSpD)

rowData(spe)$xSpD_enrich <- NULL
rowData(spe)$xSpD_enrich <- enrichment_top$xSpD[match(rownames(spe), enrichment_top$gene)] 
table(rowData(spe)$xSpD_enrich)

dotplot_enrichment <- spe |>
    scDotPlot(features = enrichment_top$gene,
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "xSpD_enrich",
              scale = TRUE,
              annoColors = list("xSpD" = SpX_colors,
                                "xSpD_enrich" = SpX_colors),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_enrichment, filename = here(plot_dir, "Xenium_xSpD_dotplot_enrichment.pdf"), width = 7, height = 8)
ggsave(dotplot_enrichment, filename = here(plot_dir, "Xenium_xSpD_dotplot_enrichment.png"), width = 7, height = 8)

#### Key genes ####

erc_xSpD_key <- read_csv(here(data_dir, "xSpD_key_genes_long.csv"))

rowData(spe)$xSpD_genes <- NULL
rowData(spe)$xSpD_genes <- erc_xSpD_key$domain[match(rownames(spe), erc_xSpD_key$gene)]
table(rowData(spe)$xSpD_genes)

dotplot_key_genes <- spe |>
    scDotPlot(features = erc_xSpD_key$gene,
              group = "xSpD",
              groupAnno = "xSpD",
              featureAnno = "xSpD_genes",
              scale = TRUE,
              annoColors = list("xSpD" = SpX_colors,
                                "xSpD_genes" = c(SpX_colors, lit = "black")),
              clusterColumns = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)

ggsave(dotplot_key_genes, filename = here(plot_dir, "Xenium_xSpD_dotplot_key_genes.pdf"), width = 5, height = 7)
ggsave(dotplot_key_genes, filename = here(plot_dir, "Xenium_xSpD_dotplot_key_genes.png"), width = 5, height = 7)


# slurmjobs::job_single('23_Xenium_expression_dotplots', create_shell = TRUE, memory = '10G', command = "Rscript 23_Xenium_expression_dotplots.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
