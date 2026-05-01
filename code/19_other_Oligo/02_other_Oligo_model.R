## Louise Huuki-Myers, Dec 2025
## Subcluster Oligos in other datasets, check for Oligo 3

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("spatialLIBD")
library("scDotPlot")
library("DeconvoBuddies")
library("getopt")

data_dir <- here("processed-data", "19_other_Oligo", "02_other_Oligo_model")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "19_other_Oligo", "02_other_Oligo_model")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Import command-line parameters
scec <- matrix(
    c("dataset", "d", "1", "character", "dataset",
      "cluster", "c", "1", "character", "clustering data",
      "opc", "o", "1", "logical", "include OPCs"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt <- list()
# opt$dataset <- "spatialDLPFC"
# opt$dataset <- "spatialdACC"
# opt$cluster <- "k30"
# opt$opc = TRUE

opt$opc <- as.logical(opt$opc)

message(Sys.time(), " - Data: ",  opt$dataset, ", Cluster: ", opt$cluster, " Include OPC: ", opt$opc)

# list.files(here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster"))
    
#### Get Oligo data ####

if(opt$dataset == "spatialDLPFC"){
    # opt$ds_short <- "dlpfc"
    
    ## get DLPFC snRNA-seq data
    sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
    sce_path <- unzip(sce_path_zip, exdir = tempdir())
    sce <- HDF5Array::loadHDF5SummarizedExperiment(
        file.path(tempdir(), "sce_DLPFC_annotated")
    )
    
    table(sce$cellType_hc)
    table(sce$cellType_layer == "Oligo")
    
    sce <- sce[,!is.na(sce$cellType_layer)]
    
    if(opt$opc){
        ## subset to Oligos & OPC
        sce <- sce[,sce$cellType_layer %in% c("Oligo", "OPC")]
        
        opt$dataset <- paste0(opt$dataset, "_wOPC")
    } else {
        ## subset to just Oligos
        sce <- sce[,sce$cellType_layer == "Oligo"]
    }
    
    sce$cellType_hc <- droplevels(sce$cellType_hc)
    table(sce$cellType_hc)
    
    sce$sample_id <- sce$BrNum
    
    dataset_covars <- c("age", "sex") # including 'Position' throws error 
    
    if(opt$cluster == "o3"){
        ## original 3 Oligo clusters from cellType_hc
        sce$Oligo_anno <- paste0("dlpfc_", gsub("_0",".", sce$cellType_hc))
        
    }
} else if(opt$dataset == "spatialHPC"){
    
    library("ExperimentHub")
    
    ehub <- ExperimentHub()
    
    ## Load the HPC dataset
    myfiles <- query(ehub, "humanHippocampus2024")
    
    sce <- myfiles[["EH9606"]]
    
    table(sce$broad.cell.class)
    # Neuron  Micro  Astro  Oligo  Other 
    # 52000   3940   4298   6787   8386 
    
    ## subset to Oligos
    sce <- sce[,sce$broad.cell.class == "Oligo"]
    sce$superfine.cell.class <- droplevels(sce$superfine.cell.class)
    
    table(sce$superfine.cell.class)
    # Oligo.1 Oligo.2 
    # 3694    3093
    
    sce$sample_id <- sce$brnum
    dataset_covars <- c("age", "sex")
} else if(opt$dataset == "spatialdACC"){
   
    load("/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/processed-data/snRNA-seq/05_azimuth/sce_azimuth.Rdata", verbose = TRUE)
    #sce
    # colData(sce)
    
    # table(sce$cellType_azimuth)
    
    sce$cellType_broad <- sce$cellType_azimuth
    sce$cellType_fine <- sce$cellType_azimuth
    
    if(opt$opc){
        ## subset to Oligos & OPC
        sce <- sce[,sce$cellType_broad %in% c("Oligo", "OPC")]
        
        opt$dataset <- paste0(opt$dataset, "_wOPC")
    } else {
        ## subset to just Oligos
        sce <- sce[,sce$cellType_broad == "Oligo"]
    }
    
    sce$sample_id <- sce$brain
    sce$doubletScore <- sce$scDblFinder.score
    
    ## no covars based on https://github.com/LieberInstitute/spatialdACC/blob/591c01ec9ca448a84242f80957af363f57b058d2/code/12_spatial_registration/azimuth_spatial_registration.R#L16-L25
    dataset_covars <- NULL
    
    ## missing logcounts
    sce <- scuttle::logNormCounts(sce)
}


## add cluster annotations
if(!"Oligo_anno" %in% colnames(colData(sce))) {
    load(here("processed-data", "19_other_Oligo", "01_other_Oligo_subcluster", sprintf("walktrap_snn_%s_subclusters_%s.Rdata", opt$cluster, opt$dataset)))
    identical(sce$key, cluster_tab$key)
    sce$Oligo_anno <- cluster_tab$cluster_anno
}

## flatten OPC subclusters to just OPCs
sce$Oligo_anno <- gsub("OPC.[0-9]","OPC", sce$Oligo_anno)

table(sce$Oligo_anno)

other_oligo_colors <- create_cell_colors(cell_types = sort(unique(sce$Oligo_anno)), palette_name = "gg")

#### run cell type specific modeling ####
modeling_fn <- here(data_dir, sprintf("modeling_results_Oligo_subtype-%s_%s.rds", opt$dataset, opt$cluster))
pseudobulk_fn = here(data_dir, sprintf("sce_pseudobulk_Oligo_subtype-%s_%s.rds", opt$dataset, opt$cluster))

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

sce_pseudo <- readRDS(here(data_dir, sprintf("sce_pseudobulk_Oligo_subtype-%s_%s.rds", opt$dataset, opt$cluster)))

top_enrichment_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    sce_layer = sce_pseudo,
    gene_name = "gene_name"
)

write.csv(top_enrichment_genes, here(data_dir, sprintf("subtype_enrichment_top10_%s_%s.csv", opt$dataset, opt$cluster)), row.names = FALSE)

#### Register with ERC Oligo subclusters ####
message(Sys.time(), " - Register vs. ERC Oligo")

load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

erc_oligo_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "modeling_results_subtype-Oligo.rds"))

cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                            modeling_results = erc_oligo_modeling,
                            model_type = "enrichment",
                            top_n = 100)

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
)

## save data
cor_layer_anno <- list(cor_layer = cor_layer, anno = anno)
write_rds(cor_layer_anno, file = here(data_dir, sprintf("other_Oligo_cor_layer_anno_%s_%s.Rds", opt$dataset, opt$cluster)))

## layer stat cor plot
rownames(cor_layer)[!rownames(cor_layer) %in% names(other_oligo_colors)]

pdf(here(plot_dir, sprintf("other_Oligo_layer_stat_cor-%s-%s.pdf", opt$dataset, opt$cluster)))
layer_stat_cor_plot(cor_layer,
                    query_colors = other_oligo_colors,
                    reference_colors = Oligo_OPC_colors,
                    annotation = anno)
dev.off()

#### Register with ERC Oligo + OPC subclusters ####
message(Sys.time(), " - Register vs. ERC Oligo + OPC")
if(opt$opc){
    erc_oligoOPC_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "35.5_sn_subcluster_marker_modeling_OligoOPC", "modeling_results_subtype-OligoOPC.rds"))
    
    cor_layer_oligoOPC <- layer_stat_cor(stats = modeling_results$enrichment,
                                         modeling_results = erc_oligoOPC_modeling,
                                         model_type = "enrichment",
                                         top_n = 100)
    
    anno_oligoOPC <- annotate_registered_clusters(
        cor_stats_layer = cor_layer_oligoOPC,
        confidence_threshold = 0.5,
        cutoff_merge_ratio = 0.1
    )
    
    ## save data
    cor_layer_anno <- list(cor_layer = cor_layer_oligoOPC, anno = anno_oligoOPC)
    write_rds(cor_layer_anno, file = here(data_dir, sprintf("other_OligoOPC_cor_layer_anno_%s_%s.Rds", opt$dataset, opt$cluster)))
    
    ## layer stat cor plot
    # rownames(cor_layer_oligoOPC)[!rownames(cor_layer_oligoOPC) %in% names(other_oligo_colors)]
    all(rownames(cor_layer_oligoOPC) %in% names(other_oligo_colors))
    
    Oligo_OPC_colors2 <- Oligo_OPC_colors[!grepl("OPC", names(Oligo_OPC_colors))]
    Oligo_OPC_colors2 <- c(Oligo_OPC_colors2, OPC = "#D2B037")
    
    pdf(here(plot_dir, sprintf("other_OligoOPC_layer_stat_cor-%s-%s.pdf", opt$dataset, opt$cluster)))
    layer_stat_cor_plot(cor_layer_oligoOPC,
                        query_colors = other_oligo_colors,
                        reference_colors = Oligo_OPC_colors2,
                        annotation = anno_oligoOPC
                        )
    dev.off()
}

#### Dot plots ####
rowData(sce)$Marker <- NULL
rowData(sce)$Marker <- top_enrichment_genes$test[match(rownames(sce), top_enrichment_genes$gene)] 
table(rowData(sce)$Marker)

top_genes_plot <- unique(top_enrichment_genes$gene[top_enrichment_genes$gene %in% rownames(sce)])

all(top_genes_plot %in% rownames(sce))

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_enrichment.pdf", opt$dataset, opt$cluster)), height = 10)
sce |>
    scDotPlot(features = top_genes_plot,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "Marker",
              annoColors = list("Oligo_anno" = other_oligo_colors,
                                "Marker" = other_oligo_colors),
              scale = FALSE,
              clusterRows = FALSE,
              groupLegends = FALSE)
sce |>
    scDotPlot(features = top_genes_plot,
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

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_lit.pdf", opt$dataset, opt$cluster)))
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

# list.files(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo"))

erc_oligo_enrich <- read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_enrichment_top10_Oligo.csv"), row.names = 1)

rowData(sce)$ERC_Oligo <- NULL
rowData(sce)$ERC_Oligo <- erc_oligo_enrich$test[match(rownames(sce), erc_oligo_enrich$gene)] 
table(rowData(sce)$ERC_Oligo)

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_ERC_Oligo_enrich.pdf", opt$dataset, opt$cluster)))
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


#### key markers from ERC subtypes ####

erc_oligo_key_genes <- list(OPC = c("PDGFRA", "MEG3","OLIG2"), #OPCs
                            COP = c("GPR17"),
                            Oligo = c("MBP", "MOG", "PLP1", "CNP", "MAL"),
                            OPC_Oligo.NF = c("RBFOX1", "KCND2", "GPM6A", #OPC + Oligo.3
                                            "CNTNAP2", "NTRK3","KCNJ3"), #OPC + Oligo.3
                            Oligo.NF = c("LINGO2", "MT-CO3", "ADGRV1"),   # Oligo.3
                            Oligo.PV = c("ARHGEF3", "ADGRF5",  "CLDN5"), #Oligo.5
                            Oligo.M = c("LAMA2", "ERBB4","OPALIN", "OMG", "SEMA6D"), #Oligo.1 + Oligo.1
                            Oligo.L = c("RASGRF1","RASGRF2", "LRRC63", "ANKRD18A") #Oligo.2
)


plot_marker_express_List(sce, 
                         cellType_col = "Oligo_anno",
                         gene_list = erc_oligo_key_genes,
                         color_pal = other_oligo_colors,
                         pdf_fn = here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_Violin_Oligo_key_genes.pdf", opt$dataset, opt$cluster))
)


oligo_key_genes <- map(erc_oligo_key_genes, ~.x[.x %in% rownames(sce)])
oligo_key_genes <- AnnotationDbi::unlist2(oligo_key_genes)

rowData(sce)$oligo_key_genes <- NULL
rowData(sce)$oligo_key_genes <- names(oligo_key_genes)[match(rownames(sce), oligo_key_genes)] 
table(rowData(sce)$oligo_key_genes)

pdf(here(plot_dir, sprintf("other_Oligo_subtype_%s_%s_dotplot_key.pdf", opt$dataset, opt$cluster)))
sce |>
    scDotPlot(features = oligo_key_genes,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "oligo_key_genes",
              scale = TRUE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)

sce |>
    scDotPlot(features = oligo_key_genes,
              group = "Oligo_anno",
              groupAnno = "Oligo_anno",
              featureAnno = "oligo_key_genes",
              scale = FALSE,
              annoColors = list("Oligo_anno" = other_oligo_colors),
              clusterRows = FALSE,
              groupLegends = FALSE)
dev.off()


# slurmjobs::job_single('02_other_Oligo_model_dacc', create_shell = TRUE, memory = '25G', command = "Rscript 02_other_Oligo_model.R --dataset spatialdACC --cluster k20 --opc TRUE")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
