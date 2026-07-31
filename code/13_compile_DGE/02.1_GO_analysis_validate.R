## GO enrichment restricted to Oligo.3 DEGs with independent Xenium validation
## Companion to 02_GO_analysis.R - instead of testing all Oligo.3 DEGs, this
## restricts the "up"/"down" gene sets to genes that also validated in >=1
## Xenium context, while keeping the full tested Oligo.3 gene set (from the
## sn_fine carrier model) as the enrichment background/universe.

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("org.Hs.eg.db")
library("clusterProfiler")
library("ComplexHeatmap")

target_cluster <- "Oligo.3" # only cluster analyzed here - kept as a variable
                            # in case this script is adapted for others later

data_dir <- here("processed-data", "13_compile_DGE", "02.1_GO_analysis_validate")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "02.1_GO_analysis_validate")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load Xenium validation data ####

validation_data_types <- c("Xenium_cell_type_anno", "Xenium_cell_type_anno_SpX", "Xenium_Oligo.3_Astro")

load_DE_valid_data <- function(datatype) {
    DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s_wSN.Rds", datatype))
    stopifnot("DE validation data file not found" = file.exists(DE_data_fn))
    message("loading ", datatype, ": ", DE_data_fn)
    readRDS(DE_data_fn)
}

data_type_lookup <- tribble(
    ~data_type,                    ~data_type_short,
    "Xenium_cell_type_anno",       "cell_type",
    "Xenium_cell_type_anno_SpX",   "SpX",
    "Xenium_Oligo.3_Astro",        "Oligo.3_Nbr"
)

DE_valid_data <- validation_data_types |>
    map(load_DE_valid_data) |>
    list_rbind() |>
    mutate(cluster = ifelse(data_type == "Xenium_cell_type_anno_SpX", cluster_SpX, as.character(cluster))) |>
    left_join(data_type_lookup)

message(nrow(DE_valid_data), " rows of Xenium validation data across ", n_distinct(DE_valid_data$data_type_short), " validation contexts")

## Restrict to the target cluster using cell_type_anno, not `cluster` -
## `cluster` holds neighbor-context labels for Xenium_Oligo.3_Astro (e.g.
## "APOE_high_nnA_Astro.2", "nn_Astro.1") with no "Oligo.3" in the name, while
## cell_type_anno stays "Oligo.3" consistently across all three data types
DE_valid_target <- DE_valid_data |> filter(cell_type_anno == target_cluster, vlmf_sn_adj.P.Val < 0.05)
stopifnot("No rows found in validation data for target_cluster" = nrow(DE_valid_target) > 0)

DE_valid_target |> count(data_type_short)

## A gene counts as validated if it validates in *any* spatial domain or
## astrocyte-neighbor context - one confirming context is sufficient evidence,
## since each context is a distinct, independent test of the same DEG
validated_genes <- DE_valid_target |>
    group_by(gene_name) |>
    summarize(validated = any(validate, na.rm = TRUE), .groups = "drop")

message(
    sum(validated_genes$validated), " / ", nrow(validated_genes),
    " tested ", target_cluster, " genes validated in >=1 Xenium context"
)

# 31 / 63 tested Oligo.3 genes validated in >=1 Xenium context

#### Load sn_fine DE data: defines the Oligo.3 gene universe ####
## Background is the full tested gene set for this cluster in the sn_fine
## carrier model - NOT just the subset represented on the Xenium panel, which
## would bias enrichment toward whatever happened to be chosen for a targeted
## panel rather than reflecting the true background of expressed genes.

DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds")
stopifnot("sn_fine DE data file not found" = file.exists(DE_data_fn))
DE_data <- readRDS(DE_data_fn)

DE_data_target <- DE_data |> filter(cluster == target_cluster)
message(nrow(DE_data_target), " genes tested for ", target_cluster, " in sn_fine carrier model (defines GO universe)")

## ENTREZID lookup
entrez_search <- bitr(DE_data_target$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
entrez_search |> count(ENSEMBL) |> count(n)

DE_entrez_target <- DE_data_target |>
    left_join(entrez_search, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") |>
    filter(!is.na(ENTREZID)) |>
    left_join(validated_genes, by = "gene_name") |>
    mutate(
        validated = replace_na(validated, FALSE),
        DE_class = case_when(
            validated & vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "validated_up",
            validated & vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "validated_down",
            TRUE ~ "None"
        ),
        ## match the "." -> "-" convention used in 02_GO_analysis.R so
        ## DE_class_cluster labels are consistent across both scripts
        DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_", DE_class)
    )

DE_entrez_target |> count(DE_class_cluster)

#### Run GO ####

universe <- unique(DE_entrez_target$ENTREZID)
length(universe)

ont_list <- c("CC", "BP", "MF")
names(ont_list) <- ont_list

go_result_validate <- map(ont_list, ~ compareCluster(ENTREZID ~ DE_class_cluster,
    data = DE_entrez_target |> filter(DE_class != "None"),
    OrgDb = org.Hs.eg.db,
    fun = enrichGO,
    universe = universe,
    ont = .x, ## ALL, CC, BP, MF
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
))

saveRDS(go_result_validate, file = here(data_dir, "GO_result_sn_fine_Oligo3_validate.rds"))

## convert to table
compare_clus_validate <- map2_dfr(go_result_validate, names(go_result_validate), ~ .x@compareClusterResult |> mutate(ONTOLOGY = .y))
compare_clus_validate |> count(DE_class_cluster, ONTOLOGY)

saveRDS(compare_clus_validate, file = here(data_dir, "GO_compare_clus_sn_fine_Oligo3_validate.rds"))
write.csv(compare_clus_validate, file = here(data_dir, "GO_results_sn_fine_Oligo3_validate.csv"), row.names = FALSE)

## quick check: are BCAS1/ERBB3/LPAR1/UGT8 (validated, not previously
## highlighted in text) driving any terms alongside the genes already named?
# compare_clus_validate |> filter(grepl("BCAS1|ERBB3|LPAR1|UGT8", geneID)) |> arrange(p.adjust)

#### dot plots ####

pdf(file = here(plot_dir, "GO_dotplot_sn_fine_Oligo3_validate.pdf"), width = 8, height = 8)
walk2(
    go_result_validate, names(go_result_validate),
    ~ print(
        dotplot(.x,
            x = "DE_class_cluster",
            showCategory = 10,
            label_format = 60
        ) +
            ggtitle(paste("GO Enrichment (Xenium-validated", target_cluster, "DEGs):", .y)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
)
dev.off()

pdf(file = here(plot_dir, "GO_dotplot_sn_fine_Oligo3_validate_small.pdf"), width = 7, height = 7)
walk2(
    go_result_validate, names(go_result_validate),
    ~ print(
        dotplot(.x,
            x = "DE_class_cluster",
            showCategory = 5,
            label_format = 30
        ) +
            ggtitle(paste("GO Enrichment\nXenium-validated", target_cluster, "DEGs:", .y)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    )
)
dev.off()

## same plot restricted to terms with >=2 validated genes, to filter out
## single-gene-driven terms (the same caveat that applied to LRRK2's broad
## GO-term footprint in the full, non-validated Oligo.3 gene set)
pdf(file = here(plot_dir, "GO_dotplot_sn_fine_Oligo3_validate_2plus.pdf"), width = 8, height = 8)
walk2(go_result_validate, names(go_result_validate), function(gr, ont) {
    gr@compareClusterResult <- gr@compareClusterResult |> filter(Count >= 2)
    if (nrow(gr@compareClusterResult) == 0) {
        return(NULL)
    }
    print(
        dotplot(gr,
            x = "DE_class_cluster",
            showCategory = 10,
            label_format = 60
        ) +
            ggtitle(paste("GO Enrichment\nXenium-validated", target_cluster, "DEGs, 2+ genes:", ont)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    )
})
dev.off()


#### 

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
