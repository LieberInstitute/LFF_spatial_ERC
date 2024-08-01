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

qc_file = here("processed-data", "02_build_spe", "05_QC_spe", "spe_qc_anno_clean.csv")
if(file.exists(qc_file)){
    qc_anno_clean <- read.csv(qc_file)
    pd <- as.data.frame(colData(spe))
    
    pd <- pd |>
        left_join(qc_anno_clean |> select(qc_anno, key = spot_name)) |>
        mutate(qc_anno = factor(ifelse(is.na(qc_anno) & in_tissue, "None", qc_anno)),
               scran_qc_anno = factor(ifelse(as.logical(scran_discard), paste0("scran-", qc_anno), as.character(qc_anno)))
        )
    
    spe$qc_anno <- pd$qc_anno
    spe$scran_qc_anno <- pd$scran_qc_anno
}

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)

spe_discrete_vars = c(
    "ManualAnnotation",
    "in_tissue",
    vars[grep("^10x_", vars)],
    vars[grep("^scran_", vars)],
    "edge_spot",
    "scran_low_lib_size_edge",
    "qc_anno",
    "scran_qc_anno",
    "ss_qc_anno",
    "qc_anno_all",
    "local_outliers"
)


spe_discrete_vars <- spe_discrete_vars[spe_discrete_vars %in% vars]

spatialLIBD::run_app(
    spe = spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    spe_discrete_vars = spe_discrete_vars,
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
