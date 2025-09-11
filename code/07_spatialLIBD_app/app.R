library("spatialLIBD")
library("markdown")

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

#### load data ####
## main spe object
spe <- readRDS("spe_ERC_app.rds")

# Add gene search
rowData(spe)$gene_search <- paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id)

## k09 pseudobulk +modeling data
spe_pb <- readRDS("spe_pseudobulk-SpD.rds")

## define spatialLIBD column
spe_pb$spatialLIBD <- spe_pb$SpD

# here::here("processed-data","05_spe_correct_cluster","20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds")

modeling_results <- readRDS("modeling-results-SpD.rds")
sig_genes <- readRDS("sig_genes_SpD.rds")

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)

## Colors 
# metadata(spe)$SpD_colors

spe_discrete_vars = c(
    "SpD",
    "ManualAnnotation",
    vars[grep("^10x_", vars)],
    vars[grep("^scran_", vars)],
    "edge_spot",
    "scran_low_lib_size_edge",
    "qc_anno",
    "scran_qc_anno",
    "ss_qc_anno",
    "qc_anno_all",
    "local_outliers",
    sprintf("BayesSpace_SVGm_k%02d", 2:28),
    sprintf("BayesSpace_Markers_k%02d", c(2,11))
)

spe_discrete_vars <- spe_discrete_vars[spe_discrete_vars %in% vars]

spatialLIBD::run_app(
    spe = spe,
    sce_layer = spe_pb,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "LFF ERC, Visium, Sp09",
    spe_discrete_vars = spe_discrete_vars,
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio",
        "edge_distance"
    ),
    default_cluster = "SpD",
    docs_path = "www"
)
