library("spatialLIBD")
library("markdown")

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

#### load data ####
spe <- readRDS("spe_ERC_app.rds")
## k09 data
spe_pb_k09 <- readRDS("spe_pseudobulk-BayesSpace_PCA_Harmony_k09.rds")
modeling_results_k09 <- readRDS("modeling_results-BayesSpace_PCA_Harmony_k09.rds")

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
    "local_outliers",
    sprintf("BayesSpace_PCA_Harmony_k%02d", 2:28)
)


spe_discrete_vars <- spe_discrete_vars[spe_discrete_vars %in% vars]

spatialLIBD::run_app(
    spe = spe,
    sce_layer = spe_pb_k09,
    modeling_results = modeling_results_k09,
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
