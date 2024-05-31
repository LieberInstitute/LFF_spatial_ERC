library("spatialLIBD")
library("markdown")
library("here")

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

spe <- readRDS(here("processed-data", "02_build_spe","spe_raw.rds"))


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
        "scran_low_lib_size_edge"
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
