## Louise Huuki-Myers, April 2026 
## Run GO analysis on interaction data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")
library("rrvgo")
library("ComplexHeatmap")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type",
      "interaction", "i", "1", "character", "interaction"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$interaction = "Age"

# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"


data_dir <- here("processed-data", "13_compile_DGE", "15_GO_analysis_interaction", paste0("GO_", opt$interaction))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "15_GO_analysis_interaction", paste0("GO_", opt$interaction))
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

## data type details
if(opt$datatype == "sn_broad"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
}else if(opt$datatype == "sn_fine"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    
    cluster_colors <- cell_type_colors$anno
    cluster_levels <- names(cell_type_colors$anno)
    
    broad_cell_types <- names(cell_type_colors$broad)
    broad_cell_types <- broad_cell_types[broad_cell_types != "Other"]
    
}else if(opt$datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
}

cluster_levels <- cluster_levels[cluster_levels != "Other"]

## load DE data ##
# DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
# file.exists(DE_data_fn)
# 
# DE_data <- readRDS(DE_data_fn)

## load DE interaction data ##
DE_data_interaction_fn <- here("processed-data", "13_compile_DGE", "05_compile_DGE_interaction", opt$datatype, sprintf("DGE_results_interaction_%s_%s.Rds", opt$interaction, opt$datatype))
file.exists(DE_data_interaction_fn)

DE_data_interaction <- readRDS(DE_data_interaction_fn)

DE_data_interaction |> count(cluster)

head(DE_data_interaction)
DE_data_interaction |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster)

## ENTREZID look up
entrez_search <- bitr(unique(DE_data_interaction$gene_id), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez_interaction <- DE_data_interaction |> 
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
                                vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class )) 

DE_entrez_interaction |> count(DE_class_cluster)


DE_entrez_interaction_summary <- DE_entrez_interaction |> 
    group_by(ENTREZID, gene_name, cluster) |> 
    filter(!grepl("None", DE_class)) |> 
    summarise(DE_classes = paste(sort(DE_class_cluster), collapse = ", "),
              n = n())
    
DE_entrez_interaction |> filter(!grepl("None", DE_class)) |> count(DE_class_cluster) |> arrange(-n)
DE_entrez_interaction |> filter(!grepl("None", DE_class))  |> arrange(DE_class_cluster) |> select(gene_name, DE_class_cluster, cluster)

#### Run GO ####
universe <- unique(DE_entrez_interaction$ENTREZID)
message("Universe n genes: ", length(universe))

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

go_result <- map(ont_list, ~compareCluster(ENTREZID ~ DE_classes,
                                           # data = DE_entrez_interaction |> filter(!grepl("None", DE_class_cluster)), 
                                           data = DE_entrez_interaction_summary, 
                                           OrgDb = org.Hs.eg.db,
                                           fun = enrichGO,
                                           universe = universe,
                                           ont = .x, ##ALL,CC,BP,MF
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 0.05,
                                           # qvalueCutoff = 0.05,
                                           minGSSize = 5, # i dont think this worked...
                                           readable = TRUE))

saveRDS(go_result, file = here(data_dir, sprintf("GO_result_interaction_%s_%s.rds", opt$interaction, opt$datatype)))

# go_result <- readRDS(here(data_dir, sprintf("GO_result_%s_%s.rds", opt$interaction, opt$datatype)))

# convert to table
compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y)) 
compare_clus |> count(DE_classes, ONTOLOGY)

# Save 
saveRDS(compare_clus, file = here(data_dir, sprintf("GO_compare_clus_interaction_%s_%s.rds", opt$interaction, opt$datatype)))
write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_interaction_%s_%s.csv", opt$interaction, opt$datatype)), row.names = FALSE)

#### dot plots ####
pdf(file = here(plot_dir, sprintf("GO_dotplot_interaction_%s_%s.pdf", opt$interaction, opt$datatype)), width = 10, height = 12)
walk2(go_result, names(go_result), function(gr, ont){
    
    gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_classes))
    
    if(nrow(gr@compareClusterResult) == 0) return(NULL)
    
    print(
        dotplot(gr, 
                x = "DE_classes", 
                showCategory = 4, 
                label_format = 60)  +
            ggtitle(paste("GO Enrichment:", ont)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
    
})
dev.off()

## plot terms w/ 2+ genes
compare_clus |> filter(Count >= 2) |> count(DE_classes)

pdf(file = here(plot_dir, sprintf("GO_dotplot_2plus_interaction_%s_%s.pdf", opt$interaction, opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), function(gr, ont){
    
    gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_classes), Count >= 5)
    
    if(nrow(gr@compareClusterResult) == 0) return(NULL)
    
    print(
        dotplot(gr, 
                x = "DE_classes", 
                showCategory = 5, 
                label_format = 60)  +
            ggtitle(paste("GO Enrichment:", ont, "2+ genes")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
    
})
dev.off()

## plot terms w/ 5+ genes
compare_clus |> filter(Count >= 5) |> count(DE_classes)

pdf(file = here(plot_dir, sprintf("GO_dotplot_5plus_interaction_%s_%s.pdf", opt$interaction, opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), function(gr, ont){
    
    gr@compareClusterResult <- gr@compareClusterResult |> filter(!grepl("0", DE_classes), Count >= 5)
    
    if(nrow(gr@compareClusterResult) == 0) return(NULL)
    
    print(
        dotplot(gr, 
                x = "DE_classes", 
                showCategory = 5, 
                label_format = 60)  +
            ggtitle(paste("GO Enrichment:", ont, "5+ genes")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
    
})
dev.off()

## opposite DEGs
# pdf(file = here(plot_dir, sprintf("GO_dotplot_opposite_%s_%s.pdf", opt$interaction, opt$datatype)), width = 10, height = 10)
# walk2(go_result, names(go_result), function(gr, ont){
#     
#     gr@compareClusterResult <- gr@compareClusterResult |> 
#         filter(grepl("0", DE_classes), Count >= 2) |>
#         mutate(DE_class_cluster = gsub("_0_", ": ", DE_class_cluster))
#     
#     if(nrow(gr@compareClusterResult) == 0) return(NULL)
#     
#     print(
#         dotplot(gr, 
#                 x = "DE_class_cluster", 
#                 showCategory = 5, 
#                 label_format = 60)  +
#             ggtitle(paste("GO Enrichment:", ont, "opposite, 2+ genes")) +
#             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#     )
#     
# })
# dev.off()

# #### sn fine data - by cell type ####
# if(opt$datatype == "sn_fine"){
#   
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_%s_cell_type.pdf" , opt$interaction, opt$datatype)), width = 8, height = 10)
#     ct_dot_plots <- map(broad_cell_types, function(ct){
#         go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
#             # subset
#             gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster) & !grepl("0", Cluster))
#             
#             if(nrow(gr@compareClusterResult) == 0) return(NULL)
#             
#             # dotplot
#             print(dotplot(gr,
#                           x = "DE_classes",
#                           showCategory = 5,
#                           label_format = 60)  +
#                       ggtitle(sprintf("GO Enrichment: %s - %s", ont, ct)) +
#                       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#             )
#         })
#     })
#     dev.off()
#        
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_opposite_%s_%s_cell_type.pdf" , opt$interaction, opt$datatype)), width = 8, height = 10)
#     ct_dot_plots <- map(broad_cell_types, function(ct){
#         go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
#             # subset
#             gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster) & grepl("0", Cluster))
#             
#             if(nrow(gr@compareClusterResult) == 0) return(NULL)
#             
#             # dotplot
#             print(dotplot(gr,
#                           x = "DE_classes",
#                           showCategory = 5,
#                           label_format = 60)  +
#                       ggtitle(sprintf("GO Enrichment: Opposite %s - %s", ont, ct)) +
#                       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#             )
#         })
#     })
#     dev.off()
#     
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_5_plus_%s_%s_cell_type.pdf", opt$interaction, opt$datatype)), width = 8, height = 8)
#     ct_dot_plots <- map(broad_cell_types, function(ct){
#         go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
#             # subset
#             gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster), Count >= 5)
#             
#             if(nrow(gr@compareClusterResult) == 0) return(NULL)
#             
#             # dotplot
#             print(dotplot(gr,
#                           x = "DE_classes",
#                           showCategory = 5,
#                           label_format = 60)  +
#                       ggtitle(sprintf("GO Enrichment: %s - %s, 5+ genes", ont, ct)) +
#                       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#             )
#         })
#     })
#     dev.off()
#     
#     go_result_Oligo3 <- map(go_result, function(gr){
#         # gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl("Oligo.3", Cluster) | grepl("Astro.3", Cluster))
#         gr@compareClusterResult <- gr@compareClusterResult |> 
#             filter(grepl("Oligo.3", Cluster)) |>
#             mutate(DE_classes2 = gsub(", ", "\n", DE_classes),
#                    DE_classes3 = gsub(", Oligo.3_", "+", DE_classes))
#         return(gr)
#     })
#     
#     map(go_result_Oligo3, ~.x@compareClusterResult |> count(Cluster))
#     map(go_result_Oligo3, ~.x@compareClusterResult |> count(DE_classes3))
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_%s_Oligo.3.pdf", opt$interaction, opt$datatype)), width = 8, height = 12)
#     walk2(go_result_Oligo3, names(go_result_Oligo3), function(gr, ont){
#         
#         if(nrow(gr@compareClusterResult) == 0) return(NULL)
#         
#         print(
#             dotplot(gr, 
#                     x = "DE_classes3", 
#                     showCategory = 6, 
#                     label_format = 60)  +
#                 ggtitle(paste("GO Enrichment:", ont)) +
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#         )
#     }
#     )
#     dev.off()   
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_5plus_%s_%s_Oligo.3.pdf", opt$interaction, opt$datatype)), width = 8, height = 8)
#     walk2(go_result_Oligo3, names(go_result_Oligo3), function(gr, ont){
#         
#         gr@compareClusterResult <- gr@compareClusterResult |> filter(Count >= 5)
#         
#         if(nrow(gr@compareClusterResult) == 0) return(NULL)
#         
#         print(
#             dotplot(gr, 
#                     x = "DE_classes3", 
#                     showCategory = 5, 
#                     label_format = 60)  +
#                 ggtitle(paste("GO Enrichment:", ont, " 5+ genes")) +
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#         )
#             }
#     )
#     dev.off()
#     
# }

#### rrvgo ####

compare_clus |> count(Cluster, is.na(qvalue), qvalue < 0.05, Count == 1)

reducedTerms_list <- map(ont_list, function(o){
    ## loop ontologies
    go_analysis <- go_result[[o]]@compareClusterResult
    go_analysis$Cluster <- droplevels(go_analysis$Cluster)
    ## get clusters
    clusters <- levels(go_analysis$Cluster)
    names(clusters) <- clusters
    
    message(Sys.time(), " - ", o)

    ont_reducedTerms <- map(clusters, function(clus){
        ## calc similarity scores
        go_analysis_c <- go_analysis |>
            filter(Cluster == clus,
                   !is.na(qvalue), ## Sometimes qvalue NA when only 1 gene
                   qvalue < 0.05
                   ) 
        
        if(nrow(go_analysis_c) == 0) return(NULL)
        
        message(sprintf("--%s (%i)--", clus, nrow(go_analysis_c)))
        
        simMatrix <- calculateSimMatrix(go_analysis_c$ID,
                                        orgdb="org.Hs.eg.db",
                                        ont=o,
                                        semdata = GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = o),
                                        method="Rel")
        
        ## group terms based on similarity
        scores <- setNames(-log10(go_analysis_c$qvalue), go_analysis_c$ID)
        reducedTerms <- NA
        
        try(reducedTerms <- reduceSimMatrix(simMatrix,
                                            scores,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db"))
        return(reducedTerms)
    }
    )

    return(ont_reducedTerms)
})

## save
saveRDS(reducedTerms_list, file = here(data_dir, sprintf("GO_reduced_terms_interaction_%s_%s.rds", opt$interaction, opt$datatype)))

reducedTerms_list2 <- list_transpose(reducedTerms_list)
reducedTerms_list2 <- reducedTerms_list2[order(names(reducedTerms_list2))]

## rm empty results
reducedTerms_list2 <- reducedTerms_list2[map_lgl(reducedTerms_list2, function(rt) !all(map_lgl(rt, ~all(is.null(.x)))))]

# map_depth(reducedTerms_list2, 2, ~.x |> dplyr::count(parentTerm))

##  treemap plots
pdf(here(plot_dir, sprintf("GO_treemap_interaction_%s_%s.pdf", opt$interaction, opt$datatype)))
walk2(reducedTerms_list2, names(reducedTerms_list2), function(rt, clus_name){
    
    rt <- rt[!map_lgl(rt, is.null)]
    map2(rt, names(rt), ~try(treemapPlot(.x, title = paste(clus_name, .y))))
    
})
dev.off()

#### GO log FC heatmap ####

reducedTerms_list_long <- map(reducedTerms_list, ~do.call("rbind", .x))

reducedTerms_list_long <- reducedTerms_list_long[!map_lgl(reducedTerms_list_long, is.null)]

source(here("code", "13_compile_DGE", "GO_logFC_heatmap.R"))


## check most common parent terms
map(reducedTerms_list_long, ~.x |> count(parentTerm) |> arrange(-n) |> head())

## check most common terms
map(reducedTerms_list_long, ~.x |> count(term) |> arrange(-n) |> head())

## check GO terms w/ most genes
compare_clus |> arrange(-Count) |> head()

