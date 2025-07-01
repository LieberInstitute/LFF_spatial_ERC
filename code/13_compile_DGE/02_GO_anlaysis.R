## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "Visium"

data_dir <- here("processed-data", "13_compile_DGE", "02_GO_analysis")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "02_GO_analysis")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

## load DE data ##
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn)

DE_data |> count(cluster)

head(DE_data)
DE_data |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster)

## ENTREZID look up
entrez_search <- bitr(DE_data$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez <- DE_data |> 
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
                                vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
                                TRUE ~ "None")) 
    
DE_entrez |> count(cluster)
DE_entrez |> count(cluster, DE_class)

DE_clusters <- DE_entrez |> filter(DE_class != "None") |> pull(cluster) |> unique()
names(DE_clusters) <- DE_clusters

go_result <- map(clusters, function(c){
    
    message(Sys.time(), " - GO ", c)
    DE_entrez <- DE_entrez |> filter(cluster == c)
    universe <- DE_entrez$ENTREZID
    
    go_result <- compareCluster(ENTREZID ~ DE_class,
                         data = DE_entrez |> filter(DE_class != "None"), 
                         OrgDb = org.Hs.eg.db,
                         fun = enrichGO,
                         universe = universe,
                         ont = "ALL", ##ALL,CC,BP,MF
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,#0.5,
                         qvalueCutoff = 1,
                         readable = TRUE)
    
    return(go_result)
})



save_go_dotplot <- function(go_result, ontology_type, plot_dir, filename_prefix = "filtered_model6_k12_GO", showCategory = 10) {
    # Extract underlying result data
    go_df <- go_result@compareClusterResult
    
    # Filter data frame for ontology and qvalue < 0.05
    go_subset <- subset(go_df, ONTOLOGY == ontology_type & qvalue < 0.05) ## how is qval differnt from p.adjust?
    # go_subset <- subset(go_df, ONTOLOGY == ontology_type & p.adjust < 0.05)
    
    # Skip if no significant terms
    if (nrow(go_subset) == 0) {
        message(paste("No significant terms for", ontology_type))
        return(NULL)
    }
    
    # Create a temporary compareClusterResult object with filtered data
    go_result@compareClusterResult <- go_subset
    
    # File path
    pdf_file <- here::here(plot_dir, paste0(filename_prefix, "_", ontology_type, ".pdf"))
    
    # Save to PDF
    pdf(file = pdf_file, width = 7, height = 7)
    print(
        dotplot(go_result, x = "DE_class", showCategory = showCategory, label_format = 60) +
            ggtitle(paste("GO Enrichment:", ontology_type)) +
            theme(
                plot.title = element_text(size = 12),
                axis.text.x = element_text(angle = 70, size = 12, vjust = 0.5, hjust = 0.5),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 12)
            )
    )
    dev.off()
    message(paste("Saved:", pdf_file))
}


save_go_dotplot(go_result = go_result$Astro, ontology_type = "BP", filename_prefix = "Astro")
