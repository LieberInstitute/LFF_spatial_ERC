## Louise Huuki-Myers, July 2026
## Run GO on Mediation genes

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")
library("rrvgo")
library("ComplexHeatmap")

plot_dir <- here("plots", "22_Mediation", "06.1_Mediation_GO_input")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "22_Mediation", "06.1_Mediation_GO_input")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####

## DE data for gene universe -- restricted to genes that were even
## eligible as mediation-screen candidates (significant DEGs in
## Oligo.3/Astro.1-3), matching the pool 01_mediation_screening.R drew
## mediator/outcome candidates from. This is a narrower universe than
## 02_GO_analysis.R uses (all tested genes) -- deliberately so, since the
## appropriate GO background here is "genes that could have been
## selected", not every gene on the platform.
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds")
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn)

DEG_universe <- DE_data |>
    filter(cluster %in% c("Oligo.3", "Astro.1", "Astro.2", "Astro.3"), vlmf_adj.P.Val < 0.05) |>
    pull(gene_id) |>
    unique()

length(DEG_universe)

## universe -> ENTREZID (gene_id here is ENSEMBL, matching 02_GO_analysis.R)
universe_entrez <- bitr(DEG_universe, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") |>
    pull(ENTREZID) |>
    unique()

length(universe_entrez)

mediator_outcome <- read.delim(here("processed-data", "22_Mediation", "out-erc_astro", "mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(pair = paste0(mediator, "|", outcome))


outcome_summary <- mediator_outcome |> 
    group_by(med_cl, outcome) |>
    summarise(n_med = n(),
              mediators = paste0(mediator, collapse = "|"))


#### Build flat mediation groups: one list per (med_cl, mediator), named like "Astro.3_IL6ST" ####

mediated_groups <- mediator_outcome |>
    group_by(med_cl, outcome) |>
    summarise(n_med = n(),
              mediators = paste0(mediator, collapse = "|")) |>
    {\(df) split(df$outcome, df$mediators)}()

length(mediated_groups)
map_int(mediated_groups, length)

# ABCA8|SV2B           CHRM3     CHRM3|IL6ST            FZD8      FZD8|NETO1 FZD8|NETO1|SV2B       FZD8|SV2B           IL6ST           NETO1      NETO1|SV2B           NPTXR    NPTXR|PLPPR4 
# 1             203               2              60               1              23              41               2               2              35              21               1 
# ST18            SV2B 
# 1              99 

all_groups <- mediated_groups

## all unique genes (now spanning discovery outcomes, Xenium-validated
## outcomes, AND mediator genes), then subset per group
symbol_entrez <- bitr(unique(unlist(all_groups)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

n_unmapped <- length(unique(unlist(all_groups))) - n_distinct(symbol_entrez$SYMBOL)
message(sprintf("%d/%d unique gene symbols failed to map to ENTREZID",
                n_unmapped, length(unique(unlist(all_groups)))))
symbol_entrez |> count(SYMBOL) |> count(n)  ## check for 1:many mappings

all_groups_entrez <- map(all_groups, function(genes) {
    symbol_entrez |> filter(SYMBOL %in% genes) |> pull(ENTREZID) |> unique()
})

n_dropped <- sum(map_int(all_groups_entrez, length) == 0)
if (n_dropped > 0) {
    message(sprintf("Dropping %d/%d groups with 0 mapped ENTREZID genes: %s",
                    n_dropped, length(all_groups_entrez),
                    paste(names(all_groups_entrez)[map_int(all_groups_entrez, length) == 0], collapse = ", ")))
}
all_groups_entrez <- all_groups_entrez[map_int(all_groups_entrez, length) > 0]


#### Run GO ####
## Using compareCluster()'s named-list interface directly -- this is why
## the flat, one-deep list with descriptive names matters: each name
## becomes the "Cluster" column value in the result, same column rrvgo
## section below already expects.
ont_list <- c("CC", "BP", "MF")
names(ont_list) <- ont_list

go_result <- map(ont_list, ~compareCluster(
    all_groups_entrez,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    universe = universe_entrez,
    ont = .x,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
))

saveRDS(go_result, file = here(data_dir, "GO_result_mediation_input.rds"))

## convert to table
compare_clus <- map2_dfr(go_result, names(go_result), ~.x@compareClusterResult |> mutate(ONTOLOGY = .y))
compare_clus |> count(Cluster, ONTOLOGY)

saveRDS(compare_clus, file = here(data_dir, "GO_compare_clus_mediation_input.rds"))
write.csv(compare_clus, file = here(data_dir, "GO_results_mediation_input.csv"), row.names = FALSE)

#### dot plots ####
pdf(file = here(plot_dir, "GO_dotplot_mediation_input.pdf"), width = 10, height = 10)
walk2(go_result, names(go_result),
      ~print(
          dotplot(.x, x = "Cluster", showCategory = 3, label_format = 60) +
              ggtitle(paste("GO Enrichment (mediation groups):", .y)) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
      )
)
dev.off()

pdf(file = here(plot_dir, "GO_dotplot_mediation_2plus_input.pdf"), width = 10, height = 10)
walk2(go_result, names(go_result), function(gr, ont) {
    gr@compareClusterResult <- gr@compareClusterResult |> filter(Count >= 2)
    if (nrow(gr@compareClusterResult) == 0) return(NULL)
    print(
        dotplot(gr, x = "Cluster", showCategory = 3, label_format = 60) +
            ggtitle(paste("GO Enrichment (mediation groups, 2+ genes):", ont)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
})
dev.off()

pdf(file = here(plot_dir, "GO_dotplot_mediation_2plus_small_input.pdf"), width = 5, height = 5)
walk2(go_result, names(go_result), function(gr, ont) {
    gr@compareClusterResult <- gr@compareClusterResult |> filter(Count >= 2)
    if (nrow(gr@compareClusterResult) == 0) return(NULL)
    print(
        dotplot(gr, x = "Cluster", showCategory = 2, label_format = 30) +
            ggtitle(paste("GO Enrichment\n(mediation groups, 2+ genes):", ont)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
})
dev.off()

pdf(file = here(plot_dir, "GO_dotplot_mediation_5plus_input.pdf"), width = 10, height = 10)
walk2(go_result, names(go_result), function(gr, ont) {
    gr@compareClusterResult <- gr@compareClusterResult |> filter(Count >= 5)
    if (nrow(gr@compareClusterResult) == 0) return(NULL)
    print(
        dotplot(gr, x = "Cluster", showCategory = 3, label_format = 60) +
            ggtitle(paste("GO Enrichment (mediation groups, 5+ genes):", ont)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
})
dev.off()

#### rrvgo ####
reducedTerms_list <- map(ont_list, function(o) {
    go_analysis <- go_result[[o]]@compareClusterResult
    go_analysis$Cluster <- droplevels(go_analysis$Cluster)
    clusters <- levels(go_analysis$Cluster)
    names(clusters) <- clusters
    
    message(Sys.time(), " - ", o)
    
    ont_reducedTerms <- map(clusters, function(clus) {
        go_analysis_c <- go_analysis |> filter(Cluster == clus, !is.na(qvalue))
        if (nrow(go_analysis_c) == 0) return(NULL)
        
        message(sprintf("--%s (%i)--", clus, nrow(go_analysis_c)))
        
        simMatrix <- calculateSimMatrix(go_analysis_c$ID,
                                        orgdb = "org.Hs.eg.db",
                                        ont = o,
                                        semdata = GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = o),
                                        method = "Rel")
        
        scores <- setNames(-log10(go_analysis_c$qvalue), go_analysis_c$ID)
        reducedTerms <- NA
        try(reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db"))
        
        return(reducedTerms)
    })
    
    return(ont_reducedTerms)
})

saveRDS(reducedTerms_list, file = here(data_dir, "GO_reduced_terms_mediation.rds"))

reducedTerms_list2 <- list_transpose(reducedTerms_list)
reducedTerms_list2 <- reducedTerms_list2[order(names(reducedTerms_list2))]
reducedTerms_list2 <- reducedTerms_list2[map_lgl(reducedTerms_list2, function(rt) !all(map_lgl(rt, ~all(is.null(.x)))))]

pdf(here(plot_dir, "GO_treemap_mediation.pdf"))
walk2(reducedTerms_list2, names(reducedTerms_list2), function(rt, clus_name) {
    rt <- rt[!map_lgl(rt, is.null)]
    map2(rt, names(rt), ~try(treemapPlot(.x, title = paste(clus_name, .y))))
})
dev.off()


#### GO term heat map ####

compare_clus |> filter(grepl("ENC1", geneID))

mediation_xenium_hits


compare_clus_long <- compare_clus |>
    as_tibble() |>
    select(Cluster, ID, Description, geneID, ONTOLOGY, Count) |>
    separate_longer_delim(geneID, delim = "/") |>
    extract(Cluster, into = c("med_cl", "mediator", "reg"),
            regex = "^(.*)_([^_]+)_([^_]+)$",
            remove = FALSE) |>
    rename(outcome = geneID)

compare_clus_long |>
    # filter(mediator == "L6")
    count(med_cl, mediator, reg)


compare_clus_long_xen_hits <- mediation_xenium_hits |>
    select(med_cl_test, med_cl, mediator, outcome_cl, outcome, mediated) |>
    inner_join(compare_clus_long, relationship = "many-to-many")


compare_clus_long_xen_hits |> count(outcome)

top_valid_terms <- compare_clus |>
    filter(!grepl("Xen", Cluster) & ID %in% compare_clus_long_xen_hits$ID) |>
    ungroup() |>
    group_by(ID, Description, ONTOLOGY) |>
    summarize(p.adjust = min(p.adjust),
              n = n()) |>
    ungroup() |>
    group_by(ONTOLOGY) |>
    slice_min(p.adjust, n = 5)

compare_clus_long_xen_hits_summary <- compare_clus_long_xen_hits |> group_by(Cluster, med_cl,outcome_cl, ONTOLOGY, ID, Description) |> mutate(n_valid = n(), genes = paste0(unique(outcome), collapse = "/"))
compare_clus_long_xen_hits_summary |> arrange(-n_valid)


compare_clus_long_xen_hits_top <- compare_clus_long_xen_hits |>
    filter(ID %in% top_valid_terms$ID) 


## start with tile plot
compare_clus_long_xen_hits_top |>
    filter(Description == "synapse", outcome == "NPTXR")

# mutate(mediator = sprintf("%s(%s)", mediator, med_cl),
#        outcome = sprintf("%s(%s)", outcome, outcome_cl),
#        valid = sprintf("%s -> %s", med_cl_test, outcome_cl)) |>

compare_clus_long_xen_hits_top|>
    filter(ONTOLOGY == "CC") |>
    ggplot(aes(x = outcome, y = Description, color = mediator, shape = outcome_cl,
               group = outcome_cl)) +
    geom_tile(fill = "white", color = "grey") +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    # geom_point(position = position_dodge2(width = 0.6, preserve = "single"), size = 2.5)+
    # facet_wrap(~ONTOLOGY, ncol = 1, scales = "free_y") +
    facet_wrap(~mediator, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()