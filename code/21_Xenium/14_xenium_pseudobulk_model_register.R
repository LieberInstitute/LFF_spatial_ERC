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
library("DeconvoBuddies")

# Import command-line parameters
scec <- matrix(
    c("var", "v", "append", "character", "registration variable"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
# opt$var <- "xSpD"

data_dir <- here("processed-data", "21_Xenium", "14_xenium_pseudobulk_model_register")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "14_xenium_pseudobulk_model_register")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####

if(opt$var == "xSpD"){
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

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, sprintf("xenium_modeling_results-%s.rds", opt$var)))

# modeling_results <- readRDS(here(data_dir, sprintf("xenium_modeling_results-%s.rds", opt$var)))

#### Extract Top Layer Enrichment Genes ####
(spe_pb <- readRDS(here(data_dir, sprintf("spe_xenium_pseudobulk-%s.rds", opt$var))))

top_DEGs <- sig_genes_extract(n = 10,
                              modeling_results = modeling_results,
                              model_type = "enrichment",
                              sce_layer = spe_pb,
                              gene_name = "gene_name") 


write.csv(top_DEGs, file = here(data_dir, sprintf("xenium_enrichment_modeling_%s_top100.csv", opt$var)), row.names = FALSE)
# top_DEGs <- read.csv(here(data_dir, sprintf("xenium_enrichment_modeling_%s_top100.csv", opt$var)))


if(opt$var == "xSpD"){
    ref_name <- "Visium_SpD"
    
    #### Register to Visium SpD ####
    reference_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))
    colnames(reference_modeling$enrichment) 
    
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    ref_colors <- SpD_colors
    q_colors <- metadata(spe)$SpX_colors
    
   
} else if(opt$var == "cell_type_anno") {
    #### Register to cell type enrichemnt data ####
    ref_name <- "sn_cell_type"
    
    reference_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "29_sn_subcluster_model_pseudobulk", "sce_subcluster_modeling_results-cell_type_anno.rds"))
    
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
    
    ref_colors <- cell_type_colors$anno
    q_colors <- cell_type_colors$anno
}


cor_layer <- layer_stat_cor(stats = modeling_results$enrichment,
                            modeling_results = reference_modeling,
                            model_type = "enrichment",
                            top_n = NULL)

## fix order to match levels
ref_levels <- names(ref_colors)
ref_levels <- ref_levels[ref_levels %in% colnames(cor_layer)]

q_levels <- names(q_colors)
q_levels <- q_levels[q_levels %in% rownames(cor_layer)]

cor_layer <- cor_layer[q_levels ,ref_levels]

anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.05 ## very strict annotation merge
)


pdf(here(plot_dir, sprintf("xenium_%s_v_%s_layer_stat_cor.pdf", opt$var, ref_name)),
    height = 5 + ncol(cor_layer)/10,
    width = 5 + ncol(cor_layer)/10)
print(layer_stat_cor_plot(cor_stats_layer = cor_layer,
                          cluster_rows =TRUE,
                          cluster_columns = TRUE,
                          reference_colors = ref_colors,
                          annotation = anno,
                          query_colors = q_colors,
                          column_title = ref_name,
                          row_title = paste0("Xenium_", opt$var)
))

print(layer_stat_cor_plot(cor_stats_layer = cor_layer,
                          cluster_rows =FALSE,
                          cluster_columns = FALSE,
                          reference_colors = ref_colors,
                          annotation = anno,
                          query_colors = q_colors,
                          column_title = ref_name,
                          row_title = paste0("Xenium_", opt$var)
))
dev.off()

#### Top Enrichment gene expression violin & dotplots ####

rownames(spe) <- rowData(spe)$gene_name

if(opt$var == "xSpD"){
    
    top_DEGs <- top_DEGs |> 
        mutate(test = factor(test, levels(spe_pb$xSpD)),
               duplicated = duplicated(gene),
               anno = paste0(top, ": logFC =", round(logFC, 2)),
               cellType.target = test) |>
        group_by(gene) |> 
        mutate(n = n()) |>
        slice_max(logFC) |>
        ungroup() |>
        arrange(test, top)
    
    DeconvoBuddies::plot_marker_express_ALL(sce = spe, 
                                        stats = top_DEGs,
                                        n_genes = 6,
                                        rank_col = "top",
                                        anno_col = "anno",
                                        color_pal = metadata(spe)$SpX_colors,
                                        gene_col = "gene",
                                        cellType_col = "xSpD",
                                        pdf = here(plot_dir, "xenium_SpX_topEnrich_violin_express.pdf"))
    
    
    rowData(spe)$xSpD_marker <- NULL
    rowData(spe)$xSpD_marker <- top_DEGs$test[match(rownames(spe), top_DEGs$gene)] 
    table(rowData(spe)$xSpD_marker)
    
    pdf(here(plot_dir, "xenium_SpX_topEnrich_dotplot.pdf"), height = 11)
    
    spe |>
        scDotPlot(features = top_DEGs$gene,
                  group = "xSpD",
                  groupAnno = "xSpD",
                  featureAnno = "xSpD_marker",
                  scale = TRUE,
                  annoColors = list(xSpD = metadata(spe)$SpX_colors,
                                    xSpD_marker = metadata(spe)$SpX_colors),
                  clusterRows = FALSE,
                  clusterColumns = FALSE,
                  groupLegends = FALSE)
    
    spe |>
        scDotPlot(features = top_DEGs$gene,
                  group = "xSpD",
                  groupAnno = "xSpD",
                  featureAnno = "xSpD_marker",
                  scale = TRUE,
                  annoColors = list(xSpD = metadata(spe)$SpX_colors,
                                    xSpD_marker = metadata(spe)$SpX_colors),
                  clusterRows = TRUE,
                  clusterColumns = TRUE,
                  groupLegends = FALSE)
    
    dev.off()
    
    ## just top5
    pdf(here(plot_dir, "xenium_SpX_topEnrich_dotplot_top5.pdf"))
    spe |>
        scDotPlot(features = top_DEGs |> filter(top <=5) |> pull(gene),
                  group = "xSpD",
                  groupAnno = "xSpD",
                  featureAnno = "xSpD_marker",
                  scale = TRUE,
                  annoColors = list(xSpD = metadata(spe)$SpX_colors,
                                    xSpD_marker = metadata(spe)$SpX_colors),
                  clusterRows = FALSE,
                  clusterColumns = FALSE,
                  groupLegends = FALSE)    
    
    spe |>
        scDotPlot(features = top_DEGs |> filter(top <=5) |> pull(gene),
                  group = "xSpD",
                  groupAnno = "xSpD",
                  featureAnno = "xSpD_marker",
                  scale = TRUE,
                  annoColors = list(xSpD = metadata(spe)$SpX_colors,
                                    xSpD_marker = metadata(spe)$SpX_colors),
                  clusterRows = FALSE,
                  clusterColumns = TRUE,
                  groupLegends = FALSE)    
    
    spe |>
        scDotPlot(features = top_DEGs |> filter(top <=5) |> pull(gene),
                  group = "xSpD",
                  groupAnno = "xSpD",
                  featureAnno = "xSpD_marker",
                  scale = TRUE,
                  annoColors = list(xSpD = metadata(spe)$SpX_colors,
                                    xSpD_marker = metadata(spe)$SpX_colors),
                  clusterRows = TRUE,
                  clusterColumns = TRUE,
                  groupLegends = FALSE)
    dev.off()
    
    
}else if(opt$var == "cell_type_anno"){
    
    top_DEGs <- top_DEGs |> 
        mutate(test = factor(test, ref_levels),
               duplicated = duplicated(gene),
               anno = paste0(top, ": logFC =", round(logFC, 2))) |>
        group_by(gene) |> 
        mutate(n = n()) |>
        slice_max(logFC) |>
        ungroup() |>
        arrange(test, top) |>
        mutate(cellType.target = droplevels(test))
    
    unique(top_DEGs$test)
    
    ## TODO
    DeconvoBuddies::plot_marker_express_ALL(sce = spe, 
                                        stats = top_DEGs,
                                        n_genes = 6,
                                        rank_col = "top",
                                        anno_col = "anno",
                                        color_pal = ref_colors,
                                        cellType_col = "cell_type_anno",
                                        gene_col = "gene",
                                        pdf = here(plot_dir, "xenium_cell_type_anno_topEnrich_violin_express.pdf"))
    
    
    rowData(spe)$cell_type_marker <- NULL
    rowData(spe)$cell_type_marker <- top_DEGs$test[match(rownames(spe), top_DEGs$gene)] 
    table(rowData(spe)$cell_type_marker)
    
    top2_DEGs <- top_DEGs |> group_by(test) |> slice(1:2) ## grab top 2
    
    pdf(here(plot_dir, "xenium_cell_type_anno_topEnrich_dotplot.pdf"), height = 10, width = 10)
    
    spe |>
        scDotPlot(features = top2_DEGs$gene,
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "cell_type_marker",
                  scale = TRUE,
                  annoColors = list(cell_type_anno = ref_colors,
                                    cell_type_marker = ref_colors),
                  clusterRows = FALSE,
                  clusterColumns = FALSE,
                  groupLegends = FALSE)
    
    spe |>
        scDotPlot(features = top2_DEGs$gene,
                  group = "cell_type_anno",
                  groupAnno = "cell_type_anno",
                  featureAnno = "cell_type_marker",
                  scale = TRUE,
                  annoColors = list(cell_type_anno = ref_colors,
                                    cell_type_marker = ref_colors),
                  clusterRows = TRUE,
                  clusterColumns = TRUE,
                  groupLegends = FALSE)
    
    dev.off()
}



#### Expression of key genes ####


rownames(spe_pb) <- rowData(spe_pb)$gene_name

if(opt$var == "xSpD"){

    APOE_v_SpX <- plot_gene_express(spe_pb, genes = "APOE", category = "xSpD", plot_type = 'boxplot', plot_points = TRUE, color_pal = metadata(spe_pb)$SpX_colors)
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


