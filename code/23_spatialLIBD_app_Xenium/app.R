## Louise Huuki-Myers, May 2026
## Xenium spatialLIBD App

library("spatialLIBD")
library("markdown")

## spatialLIBD uses golem.
## Golem is a framework for building production-grade shiny applications
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

#### load data ####
## main spe object
spe <- qs2::qs_read(here::here("processed-data", "21_Xenium", "16_xenium_app_prep","spe_xenium_app.qs2"))

## k09 pseudobulk +modeling data
spe_pb <- readRDS(here::here("processed-data", "21_Xenium", "14_xenium_pseudobulk_model_register", "spe_xenium_pseudobulk-SpX.rds"))
colnames(rowData(spe_pb))[1:2] <- c('gene_id', "gene_name")
rowData(spe_pb)$gene_search <- paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id)

# rownames(spe) <- rowData(spe)$gene_id


## define spatialLIBD column
spe_pb$spatialLIBD <- spe_pb$SpX

# here::here("processed-data","05_spe_correct_cluster","20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds")

modeling_results <- readRDS(here::here("processed-data", "21_Xenium", "14_xenium_pseudobulk_model_register", "xenium_modeling_results-SpX.rds"))

sig_genes <-
    spatialLIBD::sig_genes_extract_all(
        n = nrow(spe),
        modeling_results = modeling_results,
        sce_layer = spe_pb
    )

## Quickly explore the data
vars <- colnames(colData(spe))
colnames(colData(spe)) <- vars <- gsub("X10x", "10x", vars)

## Colors 
# metadata(spe)$SpD_colors

spe_discrete_vars = c(
    "SpX",
    "ManualAnnotation",
    "spot_class",
    "first_type",
    "second_type",
    "cell_type_anno",
    "cell_type_broad"
)

spe_discrete_vars <- spe_discrete_vars[spe_discrete_vars %in% vars]

spatialLIBD::run_app(
    spe = spe,
    sce_layer = spe_pb,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "LFF ERC, Xenium",
    spe_discrete_vars = spe_discrete_vars,
    spe_continuous_vars = c(
        "sum_gex",
        "detected_gex",
        "cell_area",
        "nucleus_area"
    ),
    default_cluster = "SpX",
    docs_path = "www",
    datatype = "Xenium"
)
