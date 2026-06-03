## Louise Huuki-Myers, May 2026
## Run spatial registration on 

#### Set Up ####

library("here")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("qs2")
library("spatialLIBD")
library("scDotPlot")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("var", "v", "1", "character", "registration variable"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
# opt$var <- "SpX"

data_dir <- here("processed-data", "21_Xenium", "14_xenium_pseudobulk_model_register")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "14_xenium_pseudobulk_model_register")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####

if(opt$var == "SpX"){
    ## load SpX version of spe (NOT filtered to singlets)
    spe_fn <- here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2")
} else if(opt$var == "cell_type_anno"){
    ## load cell type version of spe (filtered to singlets ONLY)
    spe_fn <- here("processed-data", "21_Xenium", "10_xenium_cell_types","spe_xenium_cell_types.qs2")
} else {
    stop("invalid opt")
}

message(Sys.time(), " - load data from: ", basename(spe_fn))
(spe <- qs_read(spe_fn))

rownames(spe) <- rowData(spe)$gene_id

message("ncells: ", ncol(spe))

## make APOE syntatic
spe$APOE_syn <- gsub("/", ".", spe$APOE)

var_reg <- opt$var
if(opt$var == "SpX"){
    
    ## add syntacticly valid version of SpX
    spe$SpX_syn <- gsub("~", "_", spe$SpX)
    table(spe$SpX_syn)
    
    var_reg <- "SpX_syn"
} 

#### Run Spatial Registration Function ####
message(Sys.time(), " - Running Spatial Registration on: ", var_reg)
stopifnot(var_reg %in% colnames(colData(spe)))

modeling_results <-registration_wrapper(
    sce = spe,
    var_registration = var_reg,
    var_sample_id = "sample_id",
    covars = c("APOE_syn", "Age", "Anc_Afr"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name",
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, sprintf("spe_xenium_pseudobulk-%s.rds", opt$var))
)

if(opt$var == "SpX") colnames(modeling_results$enrichment) <- gsub("_Sp", "~Sp", colnames(modeling_results$enrichment))

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("xenium_modeling_results-%s.rds", opt$var)))

#### Extract Top Layer Enrichment Genes ####
(spe_pb <- readRDS(here(data_dir, sprintf("spe_xenium_pseudobulk-%s.rds", opt$var))))

top_DEGs <- sig_genes_extract(n = 10,
                              modeling_results = modeling_results,
                              model_type = "enrichment",
                              sce_layer = spe_pb,
                              gene_name = "gene_name") 

if(opt$var == "SpX") top_DEGs$test <- gsub("_", "~", top_DEGs$test)

write.csv(top_DEGs, file = here(data_dir, sprintf("xenium_enrichment_modeling_%s_top100.csv", opt$var)), row.names = FALSE)


if(opt$var == "SpX"){
    ref_name <- "visium_SpD"
    
    #### Register to Visium SpD ####
    reference_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
    colnames(reference_modeling$enrichment) <- gsub("_Sp", "~Sp", colnames(reference_modeling$enrichment))
    colnames(reference_modeling$enrichment) 
    
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    ref_colors <- SpD_colors
    q_colors <- metadata(spe)$SpX_colors
   
} else if(opt$var == "cell_type_anno") {
    ref_name <- "sn_cell_type"
    
    reference_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_anno.rds"))
    
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
    
    ref_colors <- cell_type_colors$anno
    q_colors <- cell_type_colors$anno
}



cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                            modeling_results = reference_modeling,
                            model_type = "enrichment",
                            top_n = 100)

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.6,
    cutoff_merge_ratio = 0.05 ## very strict annotation merge
)


pdf(here(plot_dir, sprintf("xenium_%s_v_%s_layer_stat_cor.pdf", opt$var, ref_name)),
    height = 7 + ncol(cor_layer)/10,
    width = 7 + ncol(cor_layer)/10)
print(layer_stat_cor_plot(cor_stats_layer = cor_layer,
                          cluster_rows =TRUE,
                          cluster_columns = TRUE,
                          reference_colors = ref_colors,
                          annotation = anno,
                          query_colors = q_colors,
                          column_title = ref_name,
                          row_title = paste0("Xenium_", opt$var)
))
dev.off()

#### Expression of key genes
library(DeconvoBuddies)

rownames(spe_pb) <- rowData(spe_pb)$gene_name

if(opt$var == "SpX"){

    APOE_v_SpX <- plot_gene_express(spe_pb, genes = "APOE", category = "SpX", plot_type = 'boxplot', plot_points = TRUE, color_pal = metadata(spe_pb)$SpX_colors)
    ggsave(APOE_v_SpX, filename = here(plot_dir, "expression_APOE_v_SpX.png"))
    
} else if(opt$var == "cell_type_anno") {

}


# slurmjobs::job_single('14_xenium_pseudobulk_model_register', create_shell = TRUE, memory = '10G', command = "Rscript 14_xenium_pseudobulk_model_register.R --var cell_type_anno")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


