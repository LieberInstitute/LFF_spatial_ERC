## Louise Huuki-Myers, Dec 2025
## Subcluster Oligos in other datasets, check for Oligo 3

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("scDotPlot")
library("DeconvoBuddies")

data_dir <- here("processed-data", "19_other_Oligo", "02_other_Oligo_model")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "02_other_Oligo_model")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
# scec <- matrix(
#     c("datatype", "d", "1", "character", "Data type"),
#     ncol = 5, byrow = TRUE
# )
# opt <- getopt(scec)

dataset <- "spatialDLPFC"
cluster <- "k30"

message(Sys.time(), " - Data:",  dataset, ", Cluster: ", cluster)
    
#### Get Oligo data ####

if(dataset == "spatialDLPFC"){
    
    ## get DLPFC snRNA-seq data
    sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
    sce_path <- unzip(sce_path_zip, exdir = tempdir())
    sce <- HDF5Array::loadHDF5SummarizedExperiment(
        file.path(tempdir(), "sce_DLPFC_annotated")
    )
    
    table(sce$cellType_hc)
    table(sce$cellType_layer == "Oligo")
    
    sce <- sce[,!is.na(sce$cellType_layer)]
    
    ## subset to Oligos
    sce <- sce[,sce$cellType_layer == "Oligo"]
    
    sce$cellType_hc <- droplevels(sce$cellType_hc)
    table(sce$cellType_hc)
    
    sce$sample_id <- sce$BrNum
    
    dataset_covars <- c("age", "sex") # including 'Position' throws error 
    
    if(cluster == "o3"){
        ## original 3 Oligo clusters from cellType_hc
        sce$Oligo_anno <- paste0("dlpfc_", gsub("_0",".", sce$cellType_hc))
        
    }
}

if(cluster == "k30"){
    load(here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster", sprintf("walktrap_snn_k%02d_subclusters_%s.Rdata", 30, dataset)))
    identical(sce$key, cluster_tab$key)
    sce$Oligo_anno <- cluster_tab$cluster_anno
}

other_oligo_colors <- create_cell_colors(unique(sce$Oligo_anno))

#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("modeling_results_Oligo_subtype-%s_%s.rds", dataset, cluster))
pseudobulk_fn = here(data_dir, sprintf("sce_pseudobulk_Oligo_subtype-%s_%s.rds", dataset, cluster))

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

sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_Oligo_subtype-%s.rds", dataset)))

top_enrichment_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s.csv", dataset)), row.names = FALSE)

#### Register with ERC Oligo subclusters ####
load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

erc_oligo_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "modeling_results_subtype-Oligo.rds"))

cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                            modeling_results = erc_oligo_modeling,
                            model_type = "enrichment",
                            top_n = 100)

pdf(here(plot_dir, sprintf("other_Oligo_layer_stat_cor-%s-%s.pdf", dataset, cluster)))
layer_stat_cor_plot(cor_layer,
                    query_colors = other_oligo_colors,
                    reference_colors = Oligo_OPC_colors)
dev.off()


#### Dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_enrichment.pdf", dataset, cluster)))
sce |>
    scDotPlot(features = top_enrichment_genes$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "Marker",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "Marker" = other_oligo_colors),
              scale = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)
sce |>
    scDotPlot(features = top_enrichment_genes$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "Marker",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "Marker" = other_oligo_colors),
              scale = TRUE,
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## Oligo lit markers
lit_markers <- list(OPC = c("PDGFRA", "CSPG4", "MAG", "CNP", "A2B5"),
                    Oligo = c("PLP1", "ZFP191", "ZFP488", "ZFP536", "SOX17", "NKX6-2", "SMARCA4", "CD82", "TFR", "MAL"),
                    premyelin_Oligo = c("SOX10", "OLIG1", "OLIG2", "NKX2-2", "CD9", "ENPP6"),
                    myelinating_Oligo = c("BMP4", "ENPP4", "ASAP", "TMEM10", "MOG")
                    # disease_associated = c("SERPINA3", "C4B", "TNFRSF1A", "IL1B", "IL33", "HMOX1", "TNF", "ERK", "ERK2"), #https://doi.org/10.1038/s41593-025-01873-x
                    # AD_risk = c("APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "SORCS1")
)

lit_markers <- map(lit_markers, ~.x[.x %in% rownames(sce)])
lit_markers <- AnnotationDbi::unlist2(lit_markers)

rowData(sce)$litMarker <- NULL
rowData(sce)$litMarker <- names(lit_markers)[match(rownames(sce), lit_markers)] 
table(rowData(sce)$litMarker)

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_lit.pdf", dataset, cluster)))
sce |>
    scDotPlot(features = lit_markers,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "litMarker",
              scale = TRUE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = lit_markers,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "litMarker",
              scale = FALSE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()

## ERC Oligo markers

list.files(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo"))

erc_oligo_enrich <- read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_enrichment_top10_Oligo.csv"), row.names = 1)

rowData(sce)$ERC_Oligo <- NULL
rowData(sce)$ERC_Oligo <- erc_oligo_enrich$test[match(rownames(sce), erc_oligo_enrich$gene)] 
table(rowData(sce)$ERC_Oligo)

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_ERC_Oligo_enrich.pdf", dataset, cluster)))
sce |>
    scDotPlot(features = erc_oligo_enrich$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "ERC_Oligo",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "ERC_Oligo" = Oligo_OPC_colors),
              scale = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)
sce |>
    scDotPlot(features = erc_oligo_enrich$gene,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "ERC_Oligo",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "ERC_Oligo" = Oligo_OPC_colors),
              scale = TRUE,
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
