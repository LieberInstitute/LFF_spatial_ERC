## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

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
      "contrast", "c", "1", "character", "contrast"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$contrast = "ancestry"
# opt$contrast = "Sex"

# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"


data_dir <- here("processed-data", "13_compile_DGE", "10_GO_analysis_contrast", paste0("GO_", opt$contrast))
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "10_GO_analysis_contrast", paste0("GO_", opt$contrast))
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

## contrast details
if(opt$contrast == "ancestry"){
    dge_dir = "05_compile_DGE_ancestry"
    contrast_levels <- c("EA", "AA")
} else if(opt$contrast == "Sex"){
    dge_dir = "09_compile_DGE_Sex"
    contrast_levels <- c("F", "M")
} 

## load DE data ##
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn)

## load DE contrast data ##
DE_data_contrast_fn <- here("processed-data", "13_compile_DGE", dge_dir, opt$datatype, sprintf("DGE_results_%s_%s.Rds", opt$contrast, opt$datatype))
file.exists(DE_data_contrast_fn)

DE_data_contrast <- readRDS(DE_data_contrast_fn)

DE_data_contrast |> count(cluster, contrast)

head(DE_data_contrast)
DE_data_contrast |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster, contrast)

## ENTREZID look up
entrez_search <- bitr(unique(DE_data_contrast$gene_id), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez <- DE_data |> 
    filter(contrast == "carrier") |>
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
                                vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_" ,DE_class, "_ALL"),
           model_type = "All") 

DE_entrez |> count(DE_class_cluster)

DE_entrez_contrast <- DE_data_contrast |> 
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    mutate(DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "up",
                                vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class , "_", gsub("carrier_", "", contrast)),
           model_type = "Specifc") |>
    bind_rows(DE_entrez)

DE_entrez_contrast |> count(DE_class_cluster)


DE_entrez_contrast_summary <- DE_entrez_contrast |> 
    group_by(model_type, ENTREZID, gene_name, cluster) |> 
    filter(!grepl("None", DE_class)) |> 
    summarise(DE_classes = paste(sort(DE_class_cluster), collapse = ", "),
              n = n())

DE_entrez_contrast_summary |> arrange(-n)

DE_entrez_contrast_summary |> ungroup() |> count(model_type, cluster, DE_classes) 

if(opt$datatype == "sn_fine"){
    DE_entrez_contrast_summary |> ungroup() |> count(model_type, cluster, DE_classes)  |> filter(cluster == "Oligo.3")
}
    
DE_entrez_contrast |> count(DE_class_cluster)
DE_entrez_contrast |> count(cluster, contrast)
DE_entrez_contrast |> count(contrast)
DE_entrez_contrast |> filter(!grepl("None", DE_class)) |> count(DE_class_cluster) |> arrange(-n)
DE_entrez_contrast |> filter(!grepl("None", DE_class))  |> arrange(DE_class_cluster) |> select(gene_name, DE_class_cluster, cluster)

# some genes have multiple entrez matches
# DE_entrez_contrast |> filter(gene_name == "GOLGA6L10", cluster == "Astro") |> select(gene_name, gene_id, vlmf_adj.P.Val, DE_class_cluster, cluster, ENTREZID) |> unique()

#### Check contrast directions ####

# pval_contrast <- paste0("vlmf_P.Value_", contrast_levels)
# pval_adj_contrast <- paste0("vlmf_adj.P.Val_", contrast_levels)
# lfc_contrast <- paste0("vlmf_logFC_", contrast_levels)

# DE_entrez_contrast_op <- DE_data_contrast |>
#     select(gene_id, gene_name, cluster, contrast, vlmf_logFC, vlmf_P.Value, vlmf_adj.P.Val) |>
#     mutate(contrast = gsub("carrier_", "", contrast)) |>
#     pivot_wider(names_from = "contrast", values_from = c("vlmf_logFC", "vlmf_P.Value", "vlmf_adj.P.Val")) |>
#     mutate(DE_class = case_when(!!sym(lfc_contrast[[1]]) > 0 & ## up F or AA
#                                     !!sym(pval_contrast[[1]]) < 0.10 & 
#                                     !!sym(lfc_contrast[[2]]) < 0 & 
#                                     !!sym(pval_adj_contrast[[2]]) < 0.05 ~ sprintf("0_up%s_Down%s", contrast_levels[[1]], contrast_levels[[2]]),
#                                 !!sym(lfc_contrast[[1]]) > 0 & ## both UP
#                                     !!sym(pval_contrast[[1]]) < 0.10 & 
#                                     !!sym(lfc_contrast[[2]]) > 0 & 
#                                     !!sym(pval_adj_contrast[[2]]) < 0.05 ~ sprintf("0_up%s_Up%s", contrast_levels[[1]], contrast_levels[[2]]),
#                                 !!sym(lfc_contrast[[1]]) < 0 & ## down F or EA
#                                     !!sym(pval_contrast[[1]]) < 0.10 & 
#                                     !!sym(lfc_contrast[[2]]) > 0 & 
#                                     !!sym(pval_adj_contrast[[2]]) < 0.05 ~ sprintf("0_down%s_Up%s", contrast_levels[[1]], contrast_levels[[2]]),
#                                 !!sym(lfc_contrast[[1]]) < 0 & ## both Down
#                                     !!sym(pval_contrast[[1]]) < 0.10 & 
#                                     !!sym(lfc_contrast[[2]]) < 0 & 
#                                     !!sym(pval_adj_contrast[[2]]) < 0.05 ~ sprintf("0_down%s_Down%s", contrast_levels[[1]], contrast_levels[[2]]),
#                                 TRUE ~ "None"),
#            DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class)) |>
#     left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many")


# DE_entrez_contrast_op |> count(!!sym(pval_contrast[[1]]) < 0.10,
#                       !!sym(lfc_contrast[[1]]) < 0,
#                       !!sym(lfc_contrast[[2]]) < 0,
#                       !!sym(pval_adj_contrast[[2]]) < 0.05)

# DE_entrez_contrast_op |> count(DE_class)
# 
# DE_entrez_contrast_op |> filter(DE_class != "None") |> count(cluster, DE_class_cluster)
# 
# DE_opp_tile <- DE_entrez_contrast_op |>
#     count(cluster, DE_class) |>
#     ggplot(aes(x = DE_class, y = cluster, fill = n)) +
#     geom_tile() +
#     geom_text(aes(label = n, color = DE_class)) +
#     theme_bw()
# 
# ggsave(DE_opp_tile, filename = here(plot_dir, "DEG_opposite_contrast_tile.png"))
# 
# DE_entrez_contrast_op |> filter(DE_class != "None")

# DE_entrez_contrast |> filter(gene_name == "FOXO3") |> select(gene_name, contrast, cluster, vlmf_t, vlmf_adj.P.Val, vlmf_logFC) |> arrange(vlmf_adj.P.Val)

# DE_clusters <- DE_entrez_contrast |> filter(DE_class != "None") |> pull(cluster) |> unique()
# names(DE_clusters) <- DE_clusters

# DE_entrez_contrast2 <- DE_entrez_contrast |>
#     select(ENTREZID, gene_id, gene_name, cluster, contrast, DE_class, DE_class_cluster) |>
#     bind_rows(DE_entrez_contrast_op) |>
#     filter(DE_class != "None")
# 
# DE_entrez_contrast2 |> count(DE_class, is.na(ENTREZID))
# DE_entrez_contrast2 |> count(DE_class_cluster) |> print(n=50)

#### Run GO ####
universe <- unique(DE_entrez_contrast$ENTREZID)
message("Universe n genes: ", length(universe))

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

go_result <- map(ont_list, ~compareCluster(ENTREZID ~ DE_classes,
                                           # data = DE_entrez_contrast |> filter(!grepl("None", DE_class_cluster)), 
                                           data = DE_entrez_contrast_summary, 
                                           OrgDb = org.Hs.eg.db,
                                           fun = enrichGO,
                                           universe = universe,
                                           ont = .x, ##ALL,CC,BP,MF
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = 0.05,
                                           # qvalueCutoff = 0.05,
                                           minGSSize = 5, # i dont think this worked...
                                           readable = TRUE))

saveRDS(go_result, file = here(data_dir, sprintf("GO_result_%s_%s.rds", opt$contrast, opt$datatype)))

# go_result <- readRDS(here(data_dir, sprintf("GO_result_%s_%s.rds", opt$contrast, opt$datatype)))

# convert to table
compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y)) 
compare_clus |> count(DE_classes, ONTOLOGY)

# Save 
saveRDS(compare_clus, file = here(data_dir, sprintf("GO_compare_clus_%s_%s.rds", opt$contrast, opt$datatype)))
write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_%s_%s.csv", opt$contrast, opt$datatype)), row.names = FALSE)

## explore data

# compare_clus |> filter(DE_class_cluster == "Astro_AA_down")
# compare_clus |> filter(DE_class_cluster == "Astro_AA_up")

# compare_clus |> filter(DE_class_cluster == "Macro_AA_up", grepl("beta", Description))
# compare_clus |> filter(DE_class_cluster == "Vasc_EA_up")

#### dot plots ####
pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_%s.pdf", opt$contrast, opt$datatype)), width = 10, height = 12)
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

pdf(file = here(plot_dir, sprintf("GO_dotplot_2plus_%s_%s.pdf", opt$contrast, opt$datatype)), width = 10, height = 10)
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

pdf(file = here(plot_dir, sprintf("GO_dotplot_5plus_%s_%s.pdf", opt$contrast, opt$datatype)), width = 10, height = 10)
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
# pdf(file = here(plot_dir, sprintf("GO_dotplot_opposite_%s_%s.pdf", opt$contrast, opt$datatype)), width = 10, height = 10)
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

#### sn fine data - by cell type ####
if(opt$datatype == "sn_fine"){
  
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_%s_cell_type.pdf" , opt$contrast, opt$datatype)), width = 8, height = 10)
    ct_dot_plots <- map(broad_cell_types, function(ct){
        go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
            # subset
            gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster) & !grepl("0", Cluster))
            
            if(nrow(gr@compareClusterResult) == 0) return(NULL)
            
            # dotplot
            print(dotplot(gr,
                          x = "DE_classes",
                          showCategory = 5,
                          label_format = 60)  +
                      ggtitle(sprintf("GO Enrichment: %s - %s", ont, ct)) +
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
            )
        })
    })
    dev.off()
       
    pdf(file = here(plot_dir, sprintf("GO_dotplot_opposite_%s_%s_cell_type.pdf" , opt$contrast, opt$datatype)), width = 8, height = 10)
    ct_dot_plots <- map(broad_cell_types, function(ct){
        go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
            # subset
            gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster) & grepl("0", Cluster))
            
            if(nrow(gr@compareClusterResult) == 0) return(NULL)
            
            # dotplot
            print(dotplot(gr,
                          x = "DE_classes",
                          showCategory = 5,
                          label_format = 60)  +
                      ggtitle(sprintf("GO Enrichment: Opposite %s - %s", ont, ct)) +
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
            )
        })
    })
    dev.off()
    
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_5_plus_%s_%s_cell_type.pdf", opt$contrast, opt$datatype)), width = 8, height = 8)
    ct_dot_plots <- map(broad_cell_types, function(ct){
        go_result_ct <- map2(go_result, names(go_result), function(gr, ont){
            # subset
            gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl(ct, Cluster), Count >= 5)
            
            if(nrow(gr@compareClusterResult) == 0) return(NULL)
            
            # dotplot
            print(dotplot(gr,
                          x = "DE_classes",
                          showCategory = 5,
                          label_format = 60)  +
                      ggtitle(sprintf("GO Enrichment: %s - %s, 5+ genes", ont, ct)) +
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
            )
        })
    })
    dev.off()
    
    go_result_Oligo3 <- map(go_result, function(gr){
        # gr@compareClusterResult <- gr@compareClusterResult |> filter(grepl("Oligo.3", Cluster) | grepl("Astro.3", Cluster))
        gr@compareClusterResult <- gr@compareClusterResult |> 
            filter(grepl("Oligo.3", Cluster)) |>
            mutate(DE_classes2 = gsub(", ", "\n", DE_classes),
                   DE_classes3 = gsub(", Oligo.3_", "+", DE_classes))
        return(gr)
    })
    
    map(go_result_Oligo3, ~.x@compareClusterResult |> count(Cluster))
    map(go_result_Oligo3, ~.x@compareClusterResult |> count(DE_classes3))
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_%s_%s_Oligo.3.pdf", opt$contrast, opt$datatype)), width = 8, height = 12)
    walk2(go_result_Oligo3, names(go_result_Oligo3), function(gr, ont){
        
        if(nrow(gr@compareClusterResult) == 0) return(NULL)
        
        print(
            dotplot(gr, 
                    x = "DE_classes3", 
                    showCategory = 6, 
                    label_format = 60)  +
                ggtitle(paste("GO Enrichment:", ont)) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
        )
    }
    )
    dev.off()   
    
    pdf(file = here(plot_dir, sprintf("GO_dotplot_5plus_%s_%s_Oligo.3.pdf", opt$contrast, opt$datatype)), width = 8, height = 8)
    walk2(go_result_Oligo3, names(go_result_Oligo3), function(gr, ont){
        
        gr@compareClusterResult <- gr@compareClusterResult |> filter(Count >= 5)
        
        if(nrow(gr@compareClusterResult) == 0) return(NULL)
        
        print(
            dotplot(gr, 
                    x = "DE_classes3", 
                    showCategory = 5, 
                    label_format = 60)  +
                ggtitle(paste("GO Enrichment:", ont, " 5+ genes")) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
        )
            }
    )
    dev.off()
    
}

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
saveRDS(reducedTerms_list, file = here(data_dir, sprintf("GO_reduced_terms_%s_%s.rds", opt$contrast, opt$datatype)))

reducedTerms_list2 <- list_transpose(reducedTerms_list)
reducedTerms_list2 <- reducedTerms_list2[order(names(reducedTerms_list2))]

## rm empty results
reducedTerms_list2 <- reducedTerms_list2[map_lgl(reducedTerms_list2, function(rt) !all(map_lgl(rt, ~all(is.null(.x)))))]

# map_depth(reducedTerms_list2, 2, ~.x |> dplyr::count(parentTerm))

##  treemap plots
pdf(here(plot_dir, sprintf("GO_treemap_%s_%s.pdf", opt$contrast, opt$datatype)))
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

# go_terms_test <- go_lookup("synaptic membrane")


# parent_term_lookup("response to amyloid-beta")

go_genes <- get_go_genes(go_terms_test)

cluster_levels <- unlist(map(cluster_levels, ~paste0(.x, "_", contrast_levels)))

if(opt$contrast == "Ancestry"){

if(opt$datatype == "Visium"){
    
    go_terms_test <- go_lookup("axon hillock")
    go_terms_test <- go_lookup("signal recognition particle")
    
    go_genes <- get_go_genes(go_terms_test)
    get_go_genes("signal recognition particle")
    
    parent_term_lookup("signal recognition particle")
    parent_term_lookup("monoacylglycerol lipase activity")
    
    get_go_genes("signal recognition particle")
    get_go_genes("monoacylglycerol lipase activity")
    
    NPTX2
    
    compare_clus |> group_by(geneID) |> summarise(Description = paste(Description, collapse = ", "))
    
    # gene_name  DE_class_cluster      cluster      
    # <chr>      <chr>                 <fct>        
    
    # 9 NPTX2      WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07
    # 10 MAP2       WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07 # tau protein binding
    # 11 CAMKK1     WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07 # calmodulin binding
    # 12 ABHD2      WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07
    # 13 VSTM2L     WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07
    # 14 KCTD16     WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07
    # 15 WIF1       WM-uf~Sp09D07_AA_up   WM.uf~Sp09D07
    
    # 1 SRP9       WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    # 2 TMEM144    WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    # 3 CNTNAP4    WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    # 4 FAM102A    WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    # 5 GPR37      WM-uf~Sp09D07_AA_down WM.uf~Sp09D07 # ubiquitin ligase complex
    # 6 CDK18      WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    # 7 AC010642.2 WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    # 8 MRPL27     WM-uf~Sp09D07_AA_down WM.uf~Sp09D07
    
    # 16 GPR37      WM-uf~Sp09D07_EA_down WM.uf~Sp09D07 # ubiquitin ligase complex
    # 17 PEX14      WM~Sp09D06_EA_down    WM~Sp09D06   #peroxisome
    
    # geneID      Description                                                                                                                                       
    # <chr>       <chr>                                                                                                                                             
    # 1 ABHD2       sperm plasma membrane, acetylesterase activity, short-chain carboxylesterase activity, triacylglycerol lipase activity, monoacylglycerol lipase a…
    # 2 CAMKK1      calcium/calmodulin-dependent protein kinase activity                                                                                              
    # 3 CDK18       cyclin-dependent protein serine/threonine kinase activity, cyclin-dependent protein kinase activity                                               
    # 4 GPR37       receptor complex, ubiquitin ligase complex, positive regulation of catecholamine metabolic process, positive regulation of dopamine metabolic pro…
    # 5 KCTD16      G protein-coupled receptor activity involved in regulation of postsynaptic membrane potential, G protein-coupled neurotransmitter receptor activi…
    # 6 MAP2        axon hillock, basal dendrite, proximal dendrite, dendritic growth cone, dendrite terminus, apical dendrite, axon initial segment, dendrite cytopl…
    # 7 MAP2/CAMKK1 calmodulin binding                                                                                                                                
    # 8 PEX14       peroxisomal membrane, microbody membrane, peroxisome, microbody, fibrillar center, transporter complex, protein import into peroxisome matrix, su…
    # 9 SRP9        signal recognition particle, endoplasmic reticulum targeting, signal recognition particle, 7S RNA binding, signal recognition particle binding    
    # 10 TMEM144     carbohydrate transmembrane transporter activity                                                                                                   
    # 11 WIF1        Wnt-protein binding  
    
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s_%s.pdf", opt$datatype, opt$contrast)))
    
    GO_logfc_Heatmap(get_go_DE_stats("calmodulin binding", contrast = TRUE))
    GO_logfc_Heatmap(get_go_DE_stats("axon hillock", contrast = TRUE))
    GO_logfc_Heatmap(get_go_DE_stats("signal recognition particle", contrast = TRUE))
    
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("tau protein binding", 
                                             "axon hillock", 
                                             "calmodulin binding", 
                                             "triacylglycerol lipase activity"), contrast = TRUE)) # up
    
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("ubiquitin ligase complex",
                                             "receptor complex", 
                                             "peroxisomal membrane", 
                                             "microbody membrane",
                                             "peroxisome"), contrast = TRUE)) # down
    
    dev.off()
    
    go_lookup("cyclin")
    parent_term_heatmap(search_term = c("signal recognition particle",
                                        "cyclin-dependent protein kinase activity",
                                        "monoacylglycerol lipase activity",
                                        "nuclear steroid receptor activity"), 
                        pdf_suffix = "ParentTerms")
    
    
} else if(opt$datatype == "sn_broad"){
    
    ## select GO terms
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s_%s.pdf", opt$datatype, opt$contrast)))
    
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "synaptic membrane", contrast = TRUE))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "response to amyloid-beta", contrast = TRUE))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "memory", contrast = TRUE))
    
    # inhib
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("synaptic membrane","memory"), contrast = TRUE))
    # GO_logfc_Heatmap(get_go_DE_stats_multi(go_list, contrast = TRUE))
    
    dev.off()
    
    # synaptic membrane
    parent_term_heatmap(search_term = c("transport vesicle"), 
                        pdf_suffix = "test")
    
    parent_term_heatmap(search_term = c("synaptic membrane",
                                        "memory"), 
                        pdf_suffix = "ParentTerms")
    
    
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s-TOP.pdf", opt$datatype)), height = 10, width = 10)
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("synaptic membrane", 
                                             "external encapsulating structure",
                                             "regulation of trans-synaptic signaling",
                                             "extracellular structure organization",
                                             "regulation of membrane potential")),
                     title = "cellular response to calcium ion")
    dev.off()
} else if(opt$datatype == "sn_fine"){
    
    # compare_clus |> filter(grepl("MAPT", geneID))
    # compare_clus |> filter(grepl("FOS", geneID))
    
    ## select GO terms
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s.pdf", opt$datatype)), width = 12, height = 12)
    
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "neuronal cell body"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "synaptic membrane"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "metal ion transmembrane transporter activity"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "response to calcium ion"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "oligodendrocyte differentiation"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "main axon"))
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("main axon", "neuronal cell body")))
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("memory", "cognition")))
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("central nervous system myelination", "oligodendrocyte differentiation")))
    
    dev.off()
    
    # parent_term_lookup("calcium")
    
    parent_term_heatmap(search_term = c("myelin sheath",  # down in Astro.3 + Oligo.3
                                        "myelination",
                                        "cell-cell junction",
                                        "synaptic membrane",  # Oligo.3 Up
                                        "metal ion transmembrane transporter activity",
                                        "regulation of metal ion transport",
                                        "response to calcium ion",
                                        "actin filament-based movement"), 
                        pdf_suffix = "ParentTerms",
                        height = 14, width = 12)
    
}
} else if(opt$contrast == 'Sex'){
    if(opt$datatype == "sn_broad"){
        
        # compare_clus |> filter(grepl("FOS", geneID))
        
        ## select GO terms
        pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s_%s.pdf", opt$contrast, opt$datatype)), width = 12, height = 12)
        
        GO_logfc_Heatmap(get_go_DE_stats(go_term = "synaptic membrane", contrast = TRUE))
        GO_logfc_Heatmap(get_go_DE_stats(go_term = "extracellular matrix", contrast = TRUE))
        
        dev.off()
        
        # parent_term_lookup("calcium")
        
        parent_term_heatmap(search_term = c("synaptic membrane"), 
                            pdf_suffix = sprintf("%s_ParentTerms", opt$contrast),
                            height = 14, width = 12)
        
    }else if(opt$datatype == "sn_fine"){
        
        # compare_clus |> filter(grepl("MAPT", geneID))
        # compare_clus |> filter(grepl("FOS", geneID))
        
        compare_clus |> filter(grepl("chemokine", Description))
        compare_clus |> filter(grepl("eosinophil", Description))
        
        ## select GO terms
        pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s_%s.pdf", opt$contrast, opt$datatype)), width = 12, height = 12)
        
        GO_logfc_Heatmap(get_go_DE_stats(go_term = "synaptic membrane", contrast = TRUE))
        GO_logfc_Heatmap(get_go_DE_stats(go_term = "extracellular matrix", contrast = TRUE))
        
        # GO_logfc_Heatmap(get_go_DE_stats(go_term = "eosinophil migration", contrast = TRUE))
        
        dev.off()
        
        # parent_term_lookup("calcium")
        
        parent_term_heatmap(search_term = c("synaptic membrane"), 
                            pdf_suffix = sprintf("%s_ParentTerms", opt$contrast),
                            height = 14, width = 12)
        
    }
}

# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium"),
#                                  contrast = c("ancestry", "Sex")),
#                     create_shell = TRUE,
#                     name = "10_GO_analysis_contrast",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

