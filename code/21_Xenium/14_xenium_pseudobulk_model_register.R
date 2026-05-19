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
    c("brnum", "b", "1", "character", "BrNum of selected sample"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

data_dir <- here("processed-data", "21_Xenium", "14_xenium_pseudobulk_model_register")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "14_xenium_pseudobulk_model_register")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####

spe <- qs_read(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding","spe_xenium_bansky.qs2"))

## add syntacticly valid version of SpX
spe$SpX_syn <- gsub("~", "_", spe$SpX)
table(spe$SpX_syn)

## make APOE syntatic
spe$APOE <- gsub("/", "", spe$APOE)

#### Run Spatial Registration Function ####
message(Sys.time(), " - Running Spatial Registration on: SpX_syn")

modeling_results <-registration_wrapper(
    sce = spe,
    var_registration = "SpX_syn",
    var_sample_id = "sample_id",
    covars = c("APOE", "Sex", "Age", "Anc_Afr"),
    gene_ensembl = "ID",
    gene_name = "Symbol",
    min_ncells = 10,
    pseudobulk_rds_file = here(data_dir, "spe_xenium_pseudobulk-SpX.rds")
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, "xenium_modeling_results-SpX.rds"))

#### Extract Top Layer Enrichment Genes ####
(spe_pb <- readRDS(here(data_dir, "spe_xenium_pseudobulk-SpX.rds")))

top_DEGs <- sig_genes_extract(n = 10,
                              modeling_results = modeling_results,
                              model_type = "enrichment",
                              sce_layer = spe_pb,
                              gene_name = "Symbol") 

top_DEGs$test <- gsub("_", "~", top_DEGs$test)

write.csv(top_DEGs, file = here(data_dir, "enrichment_modeling_SpX_top100.csv"), row.names = FALSE)


#### Register to Visium SpD ####

layer_modeling_visium <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
colnames(layer_modeling_visium$enrichment) <- gsub("_Sp", "~Sp", colnames(layer_modeling_visium$enrichment))
colnames(layer_modeling_visium$enrichment) 

cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                            modeling_results = layer_modeling_visium,
                            model_type = "enrichment",
                            top_n = 100)


pdf(here(plot_dir, "xenium_SpX_v_visium_SpD_layer_stat_cor.pdf"))
print(layer_stat_cor_plot(cor_stats_layer = cor_layer,
                          cluster_rows =TRUE,
                          cluster_columns = TRUE))
    # ,
    # reference_colors = layer_colors[[ref]],
    # annotation = anno[[ref]],
    # query_colors = cell_type_colors$anno
dev.off()


