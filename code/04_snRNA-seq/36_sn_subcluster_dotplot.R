## Louise Huuki-Myers, July 2025
## Make dotplots of select gene sets

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("scDotPlot")
library("tidyverse")

# data_dir <- here("processed-data", "04_snRNA-seq", "36_sn_subcluster_dotplot", celltype)
# if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "36_sn_subcluster_dotplot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors

#### Sex gene dot plot ####

table(as.list(rowRanges(sce)[["seqnames"]]))
sex_check <- list(male = c("SRY", "RPS4Y1", "RPS4Y2", "DDX3Y", "KDM5D", "UTY", "ZFY", "EIF1AY", "USP9Y", "TSPY1"),
                  female = c("XIST", "TSIX", "KDM6A", "EIF2S3X", "RPS4X"))

sex_check <- map(sex_check, ~.x[.x %in% rownames(sce)])

sex_check_tb <- tibble(gene = unlist(sex_check), 
                       sex_check = unlist(map2(names(sex_check), sex_check, ~rep(.x, length(.y)))))

# Autosomal, Hormone-Responsive
# "FOXP3", "ESR1", "AR", "TSHR", "PRLR"

rowData(sce)$sex_check <- sex_check_tb$sex_check[match(rownames(sce), sex_check_tb$gene)] 
table(rowData(sce)$sex_check)


pdf(here(plot_dir, sprintf("sn_sample_dotplot_Sex_genes.pdf")))
sce |>
    scDotPlot(features = unlist(sex_check),
              group = "sample_id",
              groupAnno = "Sex",
              featureAnno = "sex_check",
              # scale = TRUE,
              annoColors = list("Sex" = sex_colors),
              clusterRows = FALSE,
              groupLegends = FALSE
              )
dev.off()

#### lit marker dot plot ####

lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    filter(gene_name %in% rowData(sce)$gene_name) |>
    # mutate(cell_type_broad = factor(cell_type_broad, levels = levels(sce$cell_type_broad))) |>
    arrange(cell_type_broad)


lit_markers |> count(cell_type_broad)
lit_markers |> count(gene_name) |> filter(n > 1)
lit_markers |> filter(gene_name == "CLDN5")


pdf(here(plot_dir, sprintf("sn_cell_type_anno_dotplot_lit_genes.pdf")))
sce |>
    scDotPlot(features = unique(lit_markers$gene_name),
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              # featureAnno = "sex_check",
              # scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = TRUE,
              groupLegends = FALSE
    )
dev.off()


#### MeanRatio dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_MeanRatio_genes$cellType.target[match(rownames(sce), top_MeanRatio_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_MeanRatio.pdf", celltype)))
sce |>
    scDotPlot(features = top_MeanRatio_genes$gene,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "Marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno,
                                "Marker" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

#### Oligo OPC marker genes ####

## from https://www.biocompare.com/Editorial-Articles/590587-A-Guide-to-Oligodendrocyte-Markers/
oligo_lit_markers <- list(OPC = c("PDGFRA", "CSPG4", "MAG", "CNP", "A2B5"),
                    Oligo = c("PLP1", "ZFP191", "ZFP488", "ZFP536", "SOX17", "NKX6-2", "SMARCA4", "CD82", "TFR", "MAL"),
                    premyelin_Oligo = c("SOX10", "OLIGO1", "OLIGO2", "NKX2-2", "CD9"),
                    myelinating_Oligo = c("BMP4", "ENPP4", "ASAP", "TMEM10", "MOG"),
                    disease_associated = c("SERPINA3", "C4B", "TNFRSF1A", "IL1B", "IL33", "HMOX1", "TNF", "ERK", "ERK2"), #https://doi.org/10.1038/s41593-025-01873-x
                    AD_risk = c("APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1")
)


oligo_lit_markers <- map(oligo_lit_markers, ~.x[.x %in% rownames(sce)])

oligo_lit_markers <- AnnotationDbi::unlist2(oligo_lit_markers)

rowData(sce)$oligo_marker <- NULL
rowData(sce)$oligo_marker <- names(oligo_lit_markers)[match(rownames(sce), oligo_lit_markers)] 
table(rowData(sce)$oligo_marker)

pdf(here(plot_dir, "sn_subtype_OligoOPC5_dotplot_lit.pdf"))
sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_anno == "OPC.5"] |>
    scDotPlot(features = oligo_lit_markers,
              group = "cell_type_anno",
              groupAnno = "cell_type_anno",
              featureAnno = "oligo_marker",
              scale = TRUE,
              annoColors = list("cell_type_anno" = cell_type_colors$anno),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

