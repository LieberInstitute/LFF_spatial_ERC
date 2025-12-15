## Louise Huuki-Myers, Dec 2025
## Subcluster Oligos in other datasets, check for Oligo 3

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("scater")
library("scran")
library("scry")

data_dir <- here("processed-data", "19_other_Oligo", "02_other_Oligo_model")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "02_other_Oligo_model")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

dataset <- "spatialDLPFC"
k <- 20

#### Get Oligo data ####

if(dataset == "spatialDLPFC"){
    
    ## get DLPFC snRNA-seq data
    sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
    sce_path <- unzip(sce_path_zip, exdir = tempdir())
    sce <- HDF5Array::loadHDF5SummarizedExperiment(
        file.path(tempdir(), "sce_DLPFC_annotated")
    )
    
    table(sce$dataset_hc)
    table(sce$dataset_layer == "Oligo")
    
    sce <- sce[,!is.na(sce$dataset_layer)]
    
    ## subset to Oligos
    sce <- sce[,sce$dataset_layer == "Oligo"]
    
    sce$dataset_hc <- droplevels(sce$dataset_hc)
    table(sce$dataset_hc)
    
    sce$sample_id <- sce$BrNum
    
    dataset_covars <- c("Position", "age", "sex")
}

load(here("processed-data", "19_other_Oligo", "02_other_Oligo_model", sprintf("walktrap_snn_k%02d_subclusters_%s.Rdata", 30, dataset)))

identical(sce$key, cluster_tab$key)

sce$Oligo_anno <- cluster_tab$cluster_anno

table(sce$Oligo_anno)

#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("modeling_results_Oligo_subtype-%s.rds", dataset))
pseudobulk_fn = here(data_dir, sprintf("sce_pseudobulk_Oligo_subtype-%s.rds", dataset))

remodel = FALSE

if(!remodel & file.exists(modeling_fn) & file.exists(pseudobulk_fn)){
    
    message(Sys.time(), " - Load modeling data")
    modeling_results <- readRDS(modeling_fn)
    
} else {
    
    message(Sys.time(), " - Running Modeling")
    
    modeling_results <-registration_wrapper(
        sce = sce,
        var_registration = "Oligo_anno",
        var_sample_id = "sample_id",
        covars = dataset_covars,
        gene_ensembl = "gene_id",
        gene_name = "gene_name",
        min_ncells = 10,
        pseudobulk_rds_file = pseudobulk_fn
    )
    
    message(Sys.time(), " - Saving Data")
    saveRDS(modeling_results, file = modeling_fn)
    
}

sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_subtype-%s.rds", dataset)))

top_enrichment_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s.csv", dataset)), row.names = FALSE)

## all signif
top100_enrichment_genes <- sig_genes_extract(
    n = 100,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

top100_enrichment_genes <- top100_enrichment_genes |>
    # filter(fdr < 0.05) |> 
    arrange(fdr) |>
    group_by(test) |>
    mutate(rank = row_number()) |>
    arrange(test)

top100_enrichment_genes |> filter(fdr < 0.05)  |> count(test)

write.csv(top100_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top100_%s.csv", dataset)), row.names = FALSE)


