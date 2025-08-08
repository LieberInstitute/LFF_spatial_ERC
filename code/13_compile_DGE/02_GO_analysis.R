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
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"

data_dir <- here("processed-data", "13_compile_DGE", "02_GO_analysis")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "02_GO_analysis")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

if(opt$datatype == "sn_broad"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
}else if(opt$datatype == "sn_fine"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$anno
    cluster_levels <- names(cell_type_colors$anno)
}else if(opt$datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
}

cluster_levels <- cluster_levels[cluster_levels != "Other"]


## load DE data ##
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn)

cluster_levels <- cluster_levels[cluster_levels %in% DE_data$cluster]

DE_data <- DE_data |> mutate(cluster = factor(cluster, levels = cluster_levels))

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
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class)) ## doesn't like .  in cluster names
    
DE_entrez |> count(cluster)
DE_entrez |> filter(DE_class != "None") |> count(DE_class_cluster)

# DE_clusters <- DE_entrez |> filter(DE_class != "None") |> pull(cluster) |> unique()
# names(DE_clusters) <- DE_clusters

#### Run GO ####
universe <- unique(DE_entrez$ENTREZID)
length(universe)

ont_list <- c("CC","BP","MF")
names(ont_list) <- ont_list

go_result <- map(ont_list, ~compareCluster(ENTREZID ~ DE_class_cluster,
                                      data = DE_entrez |> filter(DE_class != "None"), 
                                      OrgDb = org.Hs.eg.db,
                                      fun = enrichGO,
                                      universe = universe,
                                      ont = .x, ##ALL,CC,BP,MF
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      # qvalueCutoff = 0.05,
                                      readable = TRUE))

compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))
compare_clus |> count(DE_class_cluster, ONTOLOGY)

# Save data
saveRDS(compare_clus, file = here(data_dir, sprintf("GO_compare_clus_%s.rds", opt$datatype)))
write.csv(compare_clus, file = here(data_dir, sprintf("GO_results_%s.csv", opt$datatype)), row.names = FALSE)

#### dot plots ####

pdf(file = here(plot_dir, sprintf("GO_dotplot_%s.pdf", opt$datatype)), width = 10, height = 10)
walk2(go_result, names(go_result), 
      ~print(
          dotplot(.x, 
                  x = "DE_class_cluster", 
                  showCategory = 3, 
                  label_format = 60)  +
              ggtitle(paste("GO Enrichment:", .y)) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
)
dev.off()

## by cell type
# if(opt$datatype == "sn_fine"){
#     
#     pdf(file = here(plot_dir, sprintf("GO_dotplot_%s.pdf", opt$datatype)), width = 10, height = 10)
#     walk2(go_result, names(go_result), function(gr){
#         
#         ~print(
#             dotplot(.x, 
#                     x = "DE_class_cluster", 
#                     showCategory = 5, 
#                     label_format = 60)  +
#                 ggtitle(paste("GO Enrichment:", .y)) +
#                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#         )
#         
#         }
#     )
#     dev.off()
    


#### rrvgo ####

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
                   !is.na(qvalue)) ## why are some NA?
        
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


# map_depth(reducedTerms_list, 2, length)

## save
saveRDS(reducedTerms_list, file = here(data_dir, sprintf("GO_reduced_terms_%s.rds", opt$datatype)))

reducedTerms_list2 <- list_transpose(reducedTerms_list)
reducedTerms_list2 <- reducedTerms_list2[order(names(reducedTerms_list2))]

# map_depth(reducedTerms_list2, 2, length)
# map_depth(reducedTerms_list2, 2, is.null)

# ## plot treemap
# pdf(here(plot_dir, "treemap_test.pdf"))
# treemapPlot(reducedTerms, title = "Oligo")
# dev.off()

## rm empty results
reducedTerms_list2 <- reducedTerms_list2[map_lgl(reducedTerms_list2, function(rt) !all(map_lgl(rt, ~all(is.null(.x)))))]

# reducedTerms_list2[["Inhib-Vip_down"]]["CC"]

pdf(here(plot_dir, sprintf("GO_treemap_%s.pdf", opt$datatype)))
walk2(reducedTerms_list2, names(reducedTerms_list2), function(rt, clus_name){
    
    rt <- rt[!map_lgl(rt, is.null)]
    map2(rt, names(rt), ~try(treemapPlot(.x, title = paste(clus_name, .y))))
    
})
dev.off()


#### GO Heatplots ####

reducedTerms_list_long <- map(reducedTerms_list, ~do.call("rbind", .x))

## check most common parent terms
map(reducedTerms_list_long, ~.x |> count(parentTerm) |> arrange(-n) |> head())

## check most common terms
map(reducedTerms_list_long, ~.x |> count(term) |> arrange(-n) |> head())

## check GO terms w/ most genes
compare_clus |> arrange(-Count) |> head()


compare_clus |> filter(grepl("MAPT", geneID))

## find GO terms associated w/ parent terms
go_lookup <- function(search_term){
    
    hit <- any(map_lgl(reducedTerms_list_long, ~search_term %in% .x$parentTerm))
    
    if(!hit){
        stop(sprintf("GO parent term '%s' not in any reducedTerms_list", search_term))
    }
    
    go_terms <- map(reducedTerms_list_long, ~.x |> dplyr::filter(parentTerm == search_term) |> pull(go) |> unique())
    go_table <- map_dfr(go_terms, ~compare_clus |> filter(ID %in% .x))
    return(unique(go_table$Description))
    # return(go_table)
}

# go_terms_myelination <- go_lookup("myelination")
# go_terms_synaptic <- go_lookup("synaptic membrane")

parent_term_lookup <- function(search){
        map(reducedTerms_list_long, ~unique(.x$parentTerm[grep(search, .x$parentTerm)]))
}

parent_term_lookup("calcium")

## get genes matching GO terms
get_go_genes <- function(go_terms){
    
    go_genes <- compare_clus |>
        filter(Description %in% go_terms) |>
        select(ONTOLOGY, ID, Description, geneID) |>
        mutate(geneID = str_split(geneID, "/")) |>
        unnest_longer(geneID) |>
        unique()
    
    return(go_genes)
}

# go_genes <- get_go_genes(go_terms_myelination)

## gene can map parent term multiple times
# go_genes |> count(geneID) |> arrange(-n)
# go_genes |> count(Description) |> arrange(-n)

get_go_DE_stats <- function(go_term){
    
    if(!go_term %in% compare_clus$Description){
        stop(sprintf("GO term '%s' not in compare_clus", go_term))
    }
    
    go_gene <- compare_clus |>
        filter(Description == go_term) |>
        select(ONTOLOGY, ID, Description, gene_name = geneID) |>
        mutate(gene_name = str_split(gene_name, "/")) |>
        unnest_longer(gene_name) |>
        unique() 
    
    go_term_de <- go_gene |>
        left_join(DE_data |> select(cluster, gene_name, vlmf_logFC, vlmf_adj.P.Val),
                  by = "gene_name",
                  relationship = "many-to-many") |>
        group_by(cluster, gene_name) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1) |>
        mutate(signif = case_when(vlmf_adj.P.Val < 0.001 ~ "***",
                                  vlmf_adj.P.Val < 0.01 ~ "**",
                                  vlmf_adj.P.Val < 0.05 ~ "*",
                                  TRUE ~ "")
        ) 
    
    log_fc_matrix <- go_term_de |> 
        select(gene_name, cluster, vlmf_logFC) |>
        pivot_wider(names_from = "cluster", values_from = "vlmf_logFC", 
                    names_expand = TRUE) |>
        column_to_rownames("gene_name")
    
    log_fc_matrix <- log_fc_matrix[, cluster_levels]
    
    signif_matrix <- go_term_de |> 
        select(gene_name, cluster, signif) |>
        pivot_wider(names_from = "cluster", values_from = "signif", 
                    names_expand = TRUE)|>
        column_to_rownames("gene_name")
    
    signif_matrix <- signif_matrix[, cluster_levels]
    
    signif_matrix[is.na(signif_matrix)] <- ""
    
    return(list(log_fc = log_fc_matrix, 
                signif = signif_matrix,
                go_gene = go_gene))
    
}

# go_stats <- get_go_DE_stats("myelin assembly")
# go_stats <- get_go_DE_stats(go_term = "synaptic membrane")

# go_stats <- get_go_DE_stats("myelin assembly")
# go_stats <- get_go_DE_stats("myelination")
# get_go_DE_stats("protein autoprocessing")
# get_go_DE_stats("dendrite terminus")


get_go_DE_stats_multi <- function(go_list){
    
   go_stats_m <- map(go_list, get_go_DE_stats)
   go_stats_m <- map(list_transpose(go_stats_m, simplify = FALSE), ~do.call("rbind",.x))
   
   go_stats_m$go_gene <-  go_stats_m$go_gene |> 
       group_by(gene_name) |>
       summarise(Description = paste0(Description, collapse = " +\n")) |>
       arrange(Description)
   
   go_stats_m$log_fc <- go_stats_m$log_fc[go_stats_m$go_gene$gene_name,]
   go_stats_m$signif <- go_stats_m$signif[go_stats_m$go_gene$gene_name,]
   
   return(go_stats_m)
}

# go_stats_multi <- get_go_DE_stats_multi(go_list = c("myelin assembly", "protein autoprocessing", "dendrite terminus"))
# go_stats_multi <- get_go_DE_stats_multi(go_list = c("protein autoprocessing", "myelination"))

GO_logfc_Heatmap <- function(go_stats, title = NULL, cluster_rows = TRUE){
    
    logFC <- as.matrix(go_stats$log_fc)
    signif <- as.matrix(go_stats$signif)
    
    # multi_line <- any(grepl("+", go_stats$go_gene$Description))
    
    anno_table <- go_stats$go_gene |>
        # rowwise() |>
        dplyr::mutate(Description_n = ifelse(grepl("\\+", Description),
                                      Description,
                                      gsub(" ", "\n", Description))
                      ) |>
        # dplyr::mutate(Description_n = Description) |>
        column_to_rownames("gene_name")
    
    anno_table <- anno_table[rownames(logFC),]

    max_abs <- max(abs(logFC), na.rm = TRUE)
    
    my.col <- circlize::colorRamp2(
        breaks = c(-1*max_abs,0,max_abs),
        colors = c(APOE_carrier_colors[["E2+"]], "white", APOE_carrier_colors[["E4+"]])
    )
    
    
    Heatmap(as.matrix(logFC),
            col = my.col,
            name = "log(FC)",
            cluster_rows = cluster_rows,
            cluster_columns = FALSE,
            row_split = anno_table$Description,
            row_title_rot = 0,
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(signif[i, j], x, y, gp = gpar(fontsize = 10))
            },
            column_title = title)
    
}

parent_term_heatmap <- function(search_term, pdf_suffix = gsub(" ", "_", search_term), height = 10, width = 10){
    
    top_terms <- map(search_term, ~get_go_genes(go_lookup(.x))|>
                         count(Description) |>
                         arrange(-n) |>
                         head(5) |>
                         pull(Description))
    
    # TODO collapse duplicated searches
    # duplicated(top_terms)
    
    # debug
    # get_go_DE_stats_multi(top_terms[[2]])
    
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s-%s.pdf", opt$datatype, pdf_suffix)), height = height, width = width)
    
    map2(top_terms, search_term, ~print(GO_logfc_Heatmap(get_go_DE_stats_multi(.x),
                                                         title = .y)))
    dev.off()
    }


#### GO heatmap by datatype ####
if(opt$datatype == "Visium"){
    
    go_terms1 <- c("myelin assembly", 
                   "cell-cell recognition",
                   "myelination", 
                   "ensheathment of neurons", 
                   "axon ensheathment")
    
    go_terms1 %in% go_terms_myelination
    
    go_stats_multi1 <- get_go_DE_stats_multi(go_list = go_terms1)
    
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s.pdf", opt$datatype)))
    
    # GO_logfc_Heatmap(get_go_DE_stats("sperm−egg recognition"))
    GO_logfc_Heatmap(get_go_DE_stats("learning"))
    GO_logfc_Heatmap(get_go_DE_stats("myelin assembly"))
    GO_logfc_Heatmap(get_go_DE_stats("central nervous system myelination"))
    GO_logfc_Heatmap(get_go_DE_stats("oligodendrocyte differentiation"))
    GO_logfc_Heatmap(get_go_DE_stats("oligodendrocyte differentiation"))
    GO_logfc_Heatmap(go_stats_multi1)
    
    dev.off()
    
    parent_term_heatmap(search_term = c("myelination",
                                        "oligodendrocyte differentiation",
                                        "calcium ion transmembrane import into cytosol",
                                        "negative regulation of developmental growth",
                                        "dendrite terminus",
                                        "mitotic spindle midzone",
                                        "dipeptidase activity"), 
                        pdf_suffix = "ParentTerms")
    
    
} else if(opt$datatype == "sn_broad"){
    
    ## select GO terms
    pdf(here(plot_dir, sprintf("GO_logFC_heatmap_%s.pdf", opt$datatype)))

    GO_logfc_Heatmap(get_go_DE_stats(go_term = "neurotransmitter transport"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "cell fate commitment"))
    GO_logfc_Heatmap(get_go_DE_stats(go_term = "passive transmembrane transporter activity"))
    
    # inhib
    GO_logfc_Heatmap(get_go_DE_stats_multi(c("central nervous system myelination", "oligodendrocyte differentiation")))
    
    dev.off()

    parent_term_heatmap(search_term = c("synaptic membrane",
                                        "neurotransmitter transport",
                                        "cellular response to calcium ion",
                                        "learning or memory"), 
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



# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "02_GO_analysis",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

