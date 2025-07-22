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

#### Load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe

rownames(spe) <- rowData(spe)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

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


pdf(here(plot_dir, sprintf("Visium_sample_dotplot_Sex_genes.pdf")))
spe |>
    scDotPlot(features = unlist(sex_check),
              group = "sample_id",
              groupAnno = "Sex",
              featureAnno = "sex_check",
              scale = TRUE,
              annoColors = list("Sex" = sex_colors),
              clusterRows = FALSE,
              groupLegends = FALSE
              )
dev.off()


#### MeanRatio dot plots ####
rowData(spe)$Marker <- NULL
rowData(spe)$Marker <- top_MeanRatio_genes$cellType.target[match(rownames(spe), top_MeanRatio_genes$gene)] 
table(rowData(spe)$Marker)

pdf(here(plot_dir, sprintf("sn_subtype_%s_dotplot_MeanRatio.pdf", celltype)))
spe |>
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
