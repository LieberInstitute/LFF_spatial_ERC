## September 2024, Louise Huuki-Myers
## Check quality metrics and marker genes expression of Spatial Domains

library("spatialLIBD")
library("scater")
library("tidyverse")
library("DeconvoBuddies")
library("here")
library("sessioninfo")
library("readxl")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### source plotting function ####
source(here("code", "utils", "my_plot_reduced_dim.R"))

#### specify and load data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

reducedDimNames(spe)
# [1] "10x_pca"      "10x_tsne"     "10x_umap"     "PCA_p1"       "PCA_p2"       "TSNE"         "UMAP"        
# [8] "HARMONY"      "UMAP.HARMONY" "TSNE.HARMONY"

my_clusters <- sprintf("k%02d", c(9, 10, 11))
names(my_clusters) <- my_clusters

#### Load spatial registration annotations for k09, k10, k11 ####
anno_cluster_table <- map(
    my_clusters,
    ~ read_xlsx(here(
        "processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC",
        sprintf("ERC_SpD_spatial_registration_anno_summary_%s.xlsx", .x)
    ))
)

anno_cluster_table$k11 |> select(cluster, anno)

## columns to join onto colData(spe): raw cluster -> annotated SpD_k09/k10/k11
anno_cluster_table2 <- map2(
    anno_cluster_table, my_clusters,
    ~ .x |>
        mutate(SpD = fct_reorder(paste0(anno, "~", cluster), order)) |>
        select(
            !!paste0("BayesSpace_SVGm_", .y) := cluster,
            !!paste0("SpD_", .y) := SpD
        )
)

## cluster/anno/SpD lookup, used both for coloring and for annotating marker
## gene tables (which key on the raw, un-annotated cluster)
anno_cluster_table3 <- map(
    anno_cluster_table,
    ~ .x |>
        mutate(SpD = fct_reorder(paste0(anno, "~", cluster), order)) |>
        select(cluster, anno, SpD)
)

## Add SpD_k09 / SpD_k10 / SpD_k11 to colData(spe)
pd <- colData(spe) |> as.data.frame()

for (tab in anno_cluster_table2) {
    pd <- left_join(pd, tab)
}

colData(spe) <- DataFrame(pd)

table(spe$SpD_k09, spe$BayesSpace_SVGm_k09)
table(spe$SpD_k11, spe$BayesSpace_SVGm_k11)

#### Build color palette, keyed by domain name ####
SpD_colors_V4 <- c(
    "Vasc"      = "#E05AD2",  # matches vVasc
    "L1"        = "#16C72B",  # matches vL1
    "L2"        = "#021AB6",  # matches vL2
    "Inhib"     = "#C82100",  # matches vInhib
    
    ## -- L3 / LD / L5 relabeling: old-L3 -> LD, old-LD -> L5a, old-L5 -> L5b --
    "L3"        = "grey70",   # legacy name for current LD; matches vLD
    "LD"        = "grey70",   # current name; matches vLD
    "L3_Inhib"  = "#B6686F",  # k9-era compound category (LD + Inhib) - left as its own color
    "L5_LD"     = "#B1C2CE",  # ambiguous/transitional boundary category - left as its own color
    "L5a"       = "#889DF0",  # current name; matches vL5a
    "L5"        = "#0087F5",  # legacy name for current L5b; matches vL5b
    "L5b"       = "#0087F5",  # current name; matches vL5b
    
    ## -- L6 --
    "L6"        = "#40DAF2",  # matches vL6
    "L6a"       = "#83E7F7",  # k10 sublamination, no vSpD counterpart - lighter tint of vL6
    "L6b"       = "#237885",  # k10 sublamination, no vSpD counterpart - darker shade of vL6
    
    ## -- White matter --
    "WMuf"      = "#F4A460",  # matches vWMuf (unchanged)
    "WMtz"      = "#E8720C",  # legacy name for current WMim; matches vWMim (unchanged - already matched)
    "WMim"      = "#E8720C",  # vWMim
    "WMpv"      = "#E8720C",  # current name; matches vWMim
    "WM"        = "#F57A00",  # generic/undifferentiated WM catch-all, no vSpD counterpart - left as its own color
    "WMd"       = "#581009"   # matches vWMd (unchanged)
)


#' Build a color vector for one resolution's annotation table, keyed by SpD
make_SpD_colors <- function(anno_table, domain_colors = SpD_colors_V4) {
    missing_domains <- setdiff(anno_table$anno, names(domain_colors))
    if (length(missing_domains) > 0) {
        stop("No color defined for domain(s): ", paste(missing_domains, collapse = ", "))
    }

    colors_out <- domain_colors[anno_table$anno]
    names(colors_out) <- anno_table$SpD
    colors_out
}

## one color vector per resolution: SpD_colors_by_k$k09, $k10, $k11
SpD_colors_by_k <- purrr::map(anno_cluster_table3, make_SpD_colors)

#' Annotate a vector of model `test` values, which can be either single-cluster
#' tests (e.g. "Sp11D01", from anova/enrichment) or pairwise tests
#' (e.g. "Sp11D01-Sp11D02"), against a cluster/anno/SpD lookup table.
#' Pairwise tests get their annotation from each side, joined with "-"
#' (e.g. anno = "L2-L3", SpD = "L2~Sp11D01-L3~Sp11D02").
annotate_test <- function(test, anno_lookup) {
    anno_lookup <- anno_lookup |> mutate(cluster = as.character(cluster))

    tibble(test = as.character(test)) |>
        distinct() |>
        separate(test, into = c("cluster1", "cluster2"), sep = "-", fill = "right", remove = FALSE) |>
        left_join(anno_lookup |> rename(cluster1 = cluster, anno1 = anno, SpD1 = SpD), by = "cluster1") |>
        left_join(anno_lookup |> rename(cluster2 = cluster, anno2 = anno, SpD2 = SpD), by = "cluster2") |>
        mutate(
            anno = if_else(is.na(cluster2), anno1, paste0(anno1, "-", anno2)),
            SpD  = if_else(is.na(cluster2), SpD1, paste0(SpD1, "-", SpD2))
        ) |>
        select(test, anno, SpD)
}

#### Plot HARMONY reduced dims, colored by SpD_k09/k10/k11 ####
walk(c("UMAP.HARMONY", "TSNE.HARMONY"), function(dim) {
    walk(my_clusters, function(k) {
        my_plot_reduced_dim(
            spe,
            prefix = "spe",
            var_type = "cat",
            dimred = dim,
            my_var = paste0("SpD_", k),
            color_pal = SpD_colors_by_k[[k]]
        )
    })
})

#### Check QC metrics by SpD_k09/k10/k11 ####

message(Sys.time(), "- Scran Outlier Violin Plots")

map(my_clusters, function(k) {
    SpD_colors_k <- SpD_colors_by_k[[k]]
    my_var <- paste0("SpD_", k)

    pdf(here(plot_dir, sprintf("spe_erc_QC_%s.pdf", my_var)), width = 21)

    print(plotColData(spe, x = my_var, y = "sum_umi", colour_by = "scran_low_lib_size") +
        ggtitle("sum_umi") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))

    print(plotColData(spe, x = my_var, y = "sum_gene", colour_by = "scran_low_n_features") +
        ggtitle("sum_gene") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))

    print(plotColData(spe, x = my_var, y = "expr_chrM_ratio", colour_by = "scran_high_Mito_percent") +
        ggtitle("expr_chrM_ratio") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))

    dev.off()
})

pd <- as.data.frame(colData(spe))

pd_qc_long <- map(my_clusters, function(k) {
    SpD_colors_k <- SpD_colors_by_k[[k]]
    my_var <- paste0("SpD_", k)
    my_var_sym <- sym(my_var)

    pd_qc_long <- pd |>
        select(sample_id, key, sum_umi, sum_gene, expr_chrM_ratio, n_nuclei = CNmask_dark_blue, all_of(!!my_var)) |>
        mutate(n_nuclei = ifelse(n_nuclei > 100,101,n_nuclei)) |>
        pivot_longer(!c("sample_id", "key", my_var))

    SpD_qc <- pd_qc_long |>
        ggplot(aes(x = !!my_var_sym, y = value, fill = !!my_var_sym)) +
        geom_boxplot() +
        facet_wrap(~name, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = SpD_colors_k) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        )

    ggsave(SpD_qc, filename = here(plot_dir, sprintf("SpD_%s_QC.png", k)))

    return(pd_qc_long)
})

#### Check marker genes by SpD ####

rownames(spe) <- rowData(spe)$gene_name

enrich_top <- map(my_clusters, function(k) {
    spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", sprintf("spe_pseudobulk-BayesSpace_SVGm_%s.rds", k)))

    modeling_results <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", sprintf("modeling_results-BayesSpace_SVGm_%s.rds", k)))

    enrich_top <- sig_genes_extract(
        modeling_results = modeling_results,
        sce_layer = spe_pb,
        n = 10,
        model_type = "enrichment"
    )

    ## `test` holds the un-annotated cluster (matches BayesSpace_SVGm_k%02d) -
    ## join in the domain annotation as a new column
    enrich_top <- enrich_top |>
        mutate(test = as.character(test)) |>
        left_join(annotate_test(enrich_top$test, anno_cluster_table3[[k]]), by = "test")

    enrich_top_list <- map(rafalib::splitit(enrich_top$SpD), ~enrich_top$gene[.x])

    write_csv(enrich_top, file = here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", sprintf("modeling_results-BayesSpace_SVGm_%s.csv", k)))

    spe_pb$spatialLIBD <- spe_pb[[sprintf("BayesSpace_SVGm_%s", k)]]

    sig_genes_all <- sig_genes_extract_all(
        modeling_results = modeling_results,
        sce_layer = spe_pb,
        n = 10
    ) |>
        as.data.frame() |>
        mutate(test = as.character(test))

    ## `test` here includes both single-cluster tests (anova/enrichment) and
    ## pairwise tests like "Sp11D01-Sp11D02" - annotate_test() handles both,
    ## giving pairwise anno = "L2-L3" and SpD = "L2~Sp11D01-L3~Sp11D02"
    sig_genes_all <- sig_genes_all |>
        left_join(annotate_test(sig_genes_all$test, anno_cluster_table3[[k]]), by = "test")

    write_csv(sig_genes_all, file = here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", sprintf("modeling_results_ALL-BayesSpace_SVGm_%s.csv", k)))

    return(enrich_top)
})

#### compile details ####

summary_tab <- map(my_clusters, function(k) {
    my_var <- paste0("SpD_", k)

    pd_qc_median <- pd |>
        select(sum_umi, sum_gene, expr_chrM_ratio, Nmask_dark_blue, !!my_var) |>
        group_by(!!sym(my_var)) |>
        summarise(across(
            c(sum_umi, sum_gene, expr_chrM_ratio),
            list(median = median),
            na.rm = TRUE
        ))

    enrich_top_wide <- enrich_top[[k]] |>
        group_by(SpD) |>
        summarise(top_enrich_genes = paste(gene, collapse = ", ")) |>
        rename(!!my_var := SpD)

    summary_tab <- enrich_top_wide |> left_join(pd_qc_median)

    write_csv(summary_tab, file = here(data_dir, sprintf("SpD_%s_summary.csv", k)))
    return(summary_tab)
})

## marker genes from literature, plotted against annotated k09 domains
layer_marker_genes <- read.csv(here(data_dir, "layer_marker_genes.csv"))

plot_marker_express_List(
    sce = spe,
    gene_list = list(`Maynard_NatNeuro` = layer_marker_genes$gene),
    pdf_fn = here(plot_dir, "SpD_k09_layer_markers.pdf"),
    cellType_col = "SpD_k09",
    gene_name_col = "gene_name",
    color_pal = SpD_colors_by_k$k09,
    plot_points = FALSE
)

sample_order <- sort(unique(spe$BrNum))

## SVGs
nnSVG_avg <- read.csv(here("processed-data", "05_spe_correct_cluster", "13_gather_nnSVG", "nnSVG_avg.csv"))

top_svg <- nnSVG_avg$gene_name

walk(top_svg,
     ~ vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, sprintf("spe_erc-SVG_%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order
     ))

#### position ####
spe_position <- read.table(here("processed-data", "02_build_spe", "spe_position.csv"), sep = "\t", header = TRUE)

# missing order
sample_order[!sample_order %in% spe_position$sample_id]
# [1] "Br1039" "Br1556" "Br1691" "Br1706" "Br2305" "Br3974" "Br5460" "Br5529" "Br5599" "Br5832" "Br5854" "Br5941"
# [13] "Br6085" "Br6098" "Br6476"

rownames(spe_position) <- spe_position$sample_id

spe$position <- spe_position[spe$sample_id, "position"]

spe$sample_id <- paste(spe$position, spe$sample_id)

table(spe$position)

if(!dir.exists(here(plot_dir, "position"))) dir.create(here(plot_dir, "position"))

sample_order <- sort(unique(spe$sample_id))
# walk(c("MBP", "AQP4", "HPCAL1", "KRT17", "MOBP", "PCP4", "SNAP25"),
walk(c("GULP1", "COL25A1", "SATB1", "PEX5L"),
     ~ vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, "position", sprintf("spe_erc_AMY-%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order
     ))

walk(c("FN1", "NTS", "TLE1"),
     ~ vis_grid_gene(
         spe = spe,
         geneid = .x,
         pdf = here::here(plot_dir, "position", sprintf("spe_erc_HP-%s.pdf", .x)),
         assayname = "logcounts",
         point_size = 1,
         sample_order = sample_order
     ))

#### spot plots ####
if(!dir.exists(here(plot_dir, "spot_plots"))) dir.create(here(plot_dir, "spot_plots"))

sample_ids <- sort(unique(spe$sample_id))

p_list <- vis_grid_clus(
    spe = spe,
    clustervar = 'BayesSpace_SVGm_k02',
    pdf = here(plot_dir, "spot_plots", "spe_BayesSpace_SVGm_k02-ALL.pdf"),
    sort_clust = FALSE,
    spatial = FALSE,
    point_size = 1,
    sample_order = sample_ids
)


#### check k11 WM split ####

k <- "k11"

spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", sprintf("spe_pseudobulk-BayesSpace_SVGm_%s.rds", k)))
modeling_results <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", sprintf("modeling_results-BayesSpace_SVGm_%s.rds", k)))
spe_pb$spatialLIBD <- spe_pb[[sprintf("BayesSpace_SVGm_%s", k)]]

sig_genes_all <- sig_genes_extract_all(
    modeling_results = modeling_results,
    sce_layer = spe_pb,
    n = nrow(spe_pb)
) |>
    as.data.frame() |>
    mutate(test = as.character(test))

sig_genes_all <- sig_genes_all |>
    left_join(annotate_test(sig_genes_all$test, anno_cluster_table3[[k]]), by = "test")

sig_genes_all |> count(SpD) |> filter(grepl("WMpv", SpD))

# WMpv~Sp11D11-WMd~Sp11D07

# 10 Sp11D11 WMim  WMim~Sp11D11 
# 11 Sp11D07 WMd   WMd~Sp11D07

FDR_cut <- 0.01

WM_deep_pair <- sig_genes_all |>
    filter(model_type == "pairwise") |>
    filter(SpD == "WMa~Sp11D11-WMd~Sp11D07") |>  #SpD == "WMd~Sp11D07-WMim~Sp11D11"
    mutate(DE_class = case_when(fdr < FDR_cut & abs(logFC) > 1  ~ "both",
                                fdr < FDR_cut ~ paste("FDR<", FDR_cut),
                                abs(logFC) > 1 ~ "abs(logFC)>1",
                                TRUE ~ "None"))

WM_deep_pair |> filter(gene == "APLP1")

signif_colors <- c("purple", "blue", "red")
names(signif_colors) <- c("both", paste("FDR<", FDR_cut) , "abs(logFC)>1" )

library(ggrepel)

WM_deep_pair_volcano <- WM_deep_pair |>
    ggplot(aes(x = logFC, y = -log10(pval), color = DE_class)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = ifelse(DE_class != "None", gene, "")), size = 2) +
    scale_color_manual(values = signif_colors) +
    theme_bw()

ggsave(WM_deep_pair_volcano, filename = here(plot_dir, "k11_WM_deep_pair_volcano.png"))

Oligo_subtype_modeling <- readRDS(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "modeling_results_subtype-Oligo.rds"))

Oligo_subtype_top_enrich <- read_csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_enrichment_top10_Oligo.csv"))
Oligo_subtype_MeanRatio <- read_csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", "Oligo", "subtype_MeanRatio_top10_Oligo.csv"))

Oligo_subtype_check <-  Oligo_subtype_MeanRatio |> 
    select(gene, cellType.target) |>
    full_join(Oligo_subtype_top_enrich |> 
                  select(gene, cellType.target = test))

head(Oligo_subtype_modeling$enrichment)

Oligo_tstat <- Oligo_subtype_modeling$enrichment[, grepl('t_stat', colnames(Oligo_subtype_modeling$enrichment))]


WM_deep_pair_oligo <- WM_deep_pair |> left_join(Oligo_subtype_check |> select(gene, cellType.target))

WM_deep_pair_oligo |> count(cellType.target)

WM_deep_pair_volcano <- WM_deep_pair_oligo |>
    filter(!is.na(cellType.target)) |>
    ggplot(aes(x = logFC, y = -log10(pval), color = cellType.target)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = ifelse(!is.na(cellType.target), gene, "")), size = 2) +
    # scale_color_manual(values = signif_colors) +
    theme_bw()

ggsave(WM_deep_pair_volcano, filename = here(plot_dir, "k11_WM_deep_pair_volcano_Oligo_MR.png"))


## align genes by ENSEMBL ID
common_genes <- intersect(WM_deep_pair$ensembl, rownames(Oligo_tstat))
length(common_genes)  # sanity check overlap size

wm_stat   <- WM_deep_pair$stat[match(common_genes, WM_deep_pair$ensembl)]
oligo_sub <- Oligo_tstat[common_genes, , drop = FALSE]

## quick look: correlation of WM_deep_pair$stat vs each Oligo_tstat column
sapply(colnames(oligo_sub), function(col) cor(wm_stat, oligo_sub[, col]))

# t_stat_Oligo.1 t_stat_Oligo.2 t_stat_Oligo.3 t_stat_Oligo.4 t_stat_Oligo.5 
#     -0.4792996     -0.4481378      0.5423421     -0.1651901      0.3425903 


n_top <- 20

## highest and lowest WM_deep_pair genes
wm_top_high <- WM_deep_pair |> slice_max(stat, n = n_top) |> pull(ensembl)
wm_top_low  <- WM_deep_pair |> slice_min(stat, n = n_top) |> pull(ensembl)

## top 20 genes per Oligo_tstat column, keeping rownames (ensembl ids) attached
oligo_top_genes <- apply(Oligo_tstat, 2, function(col) {
    named_col <- setNames(col, rownames(Oligo_tstat))
    names(sort(named_col, decreasing = TRUE))[1:n_top]
})
# columns = Oligo_tstat columns, rows = rank 1..20 gene ids

## overlap with WM_deep_pair's top/bottom
overlap_high <- apply(oligo_top_genes, 2, function(genes) intersect(genes, wm_top_high))
overlap_low  <- apply(oligo_top_genes, 2, function(genes) intersect(genes, wm_top_low))

overlap_high
overlap_low


overlap_summary <- map_dfr(colnames(oligo_top_genes), function(col) {
    genes <- oligo_top_genes[, col]
    tibble(
        column = col,
        n_overlap_high = length(intersect(genes, wm_top_high)),
        genes_high = paste(intersect(genes, wm_top_high), collapse = ", "),
        n_overlap_low = length(intersect(genes, wm_top_low)),
        genes_low = paste(intersect(genes, wm_top_low), collapse = ", ")
    )
})

overlap_summary

# slurmjobs::job_single('11_SpD_check', create_shell = TRUE, memory = '25G', command = "Rscript 11_SpD_check.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
