library("spatialLIBD")
library("markdown")
library("tidyverse")
library("here")

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

spe <- readRDS(here("processed-data", "02_build_spe","spe_raw.rds"))

## Read in Annotations ##

qc_anno <- map_dfr(list.files(here("processed-data", "02_build_spe", "03_qc_app"), full.names = TRUE), read.csv) |>
    filter(!is.na(ManualAnnotation)) |>
    count(sample_id, spot_name, ManualAnnotation) |> 
    filter(!(spot_name == "TCTTTCCTTCGAGATA-1_Br5367" & ManualAnnotation == "out_tissue"))

qc_anno |> count(ManualAnnotation)

qc_anno |> count(sample_id, spot_name) |> filter(n>1)
# qc_anno |> filter((spot_name == "TCTTTCCTTCGAGATA-1_Br5367" & ManualAnnotation == "out_tissue"))

rownames(qc_anno) <- qc_anno$spot_name
qc_anno_col <- qc_anno[colnames(spe),]$ManualAnnotation
table(qc_anno_col)
table(is.na(qc_anno_col))
length(qc_anno_col)

spe$qc_anno <- factor(qc_anno_col)
# spe[,qc_anno$spot_name]$ManualAnnotation <- qc_anno$ManualAnnotation
# spe$ManualAnnotation <- factor(spe$ManualAnnotation)
table(spe$qc_anno)
table(spe$sample_id, spe$qc_anno)

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)
spatialLIBD::run_app(
    spe = spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    spe_discrete_vars = c(
        "ManualAnnotation",
        "in_tissue",
        vars[grep("^10x_", vars)],
        vars[grep("^scran_", vars)],
        "edge_spot",
        "scran_low_lib_size_edge",
        "qc_anno"
        # vars[grep("^SNN_k10", vars)],
        # vars[grep("^BayesSpace_harmony_", vars)]
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "edge_distance"
    ),
    default_cluster = "10x_graphclust",
    docs_path = "www"
)
