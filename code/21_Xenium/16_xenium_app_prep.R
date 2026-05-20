## Louise Huuki-Myers, May 2026
## Xenium visualizations 

#### Set Up ####

library("spatialLIBD")
library("qs2")
library("here")
library("sessioninfo")
library("tidyverse")

data_dir <- here("processed-data", "21_Xenium", "16_xenium_app_prep")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), " - Load SPE data")
spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

## Drop data we don't need for iSEE
assayNames(spe)
# [1] "counts"                      "binomial_deviance_residuals" "logcounts"
assays(spe)$counts <- NULL
assays(spe)$nucleus_normcounts <- NULL
assays(spe)$cell_normcounts <- NULL
assays(spe)$H0 <- NULL
assays(spe)$H1 <- NULL

class(logcounts(spe))
# [1] "dgCMatrix"

# sourcing official color palette
SpX_colors <- metadata(spe)$SpX_colors

## Drop metadata we don't need
metadata(spe) <- list()
spe$path <- NULL
spe$total <- NULL

reducedDims(spe)$PCA <- NULL
reducedDims(spe)$TSNE <- NULL
reducedDims(spe)$UMAP <- NULL
reducedDims(spe)$UMAP_M1_lam0.8 <- NULL

lobstr::obj_size(spe)
# 1.34 GB

colnames(rowData(spe))[1:2] <- c('gene_id', "gene_name")
rowData(spe)$gene_search <- paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id)

rownames(spe) <- rowData(spe)$gene_id

spe$key <- colnames(spe)
spe$ManualAnnotation <- spe$SpX


check_spe(spe, 
          datatype = "Xenium",
          variables = c(
              "sum_gex",  # "sum_umi",
              "detected_gex" #"sum_gene",
          ))

message(Sys.time(), " - Save qs2")
qs2::qs_save(spe, file = here(data_dir, "spe_xenium_app.qs2"))


