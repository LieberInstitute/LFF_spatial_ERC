## Louise Huuki-Myers July 2026
## All analysis summary 


#### set up ####
library("here")
library("tidyverse")
library("qs2")
library("sessioninfo")
library("MOFA2")
library("patchwork")

data_dir <- here("processed-data", "13_compile_DGE", "18_all_analysis_summary")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "18_all_analysis_summary")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
# APOE_carrier_colors
# E2+       E4+ 
#"#51B8B1" "#D97D59" 

#### Load Data ####

## load DE data ##

de_data_types <- c("sn_fine", "Visium", "Xenium_cell_type_anno", "Xenium_Oligo.3_Astro", "Xenium_SpX") 

load_DE_data <- function(datatype) {
    DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s.Rds", datatype))
    stopifnot("DE data file not found" = file.exists(DE_data_fn))
    message("loading ", datatype, ": ", DE_data_fn)
    readRDS(DE_data_fn)
}

## each Rds already carries its own data_type column, so a plain bind is fine -
## bind_rows()/list_rbind() fill NA for columns not shared across data types
## (e.g. Xenium's cell_type_anno/SpX vs sn_fine/Visium's seqnames/start/end)
DE_data <- de_data_types |>
    map(load_DE_data) |>
    list_rbind()

message(nrow(DE_data), " rows across ", n_distinct(DE_data$data_type), " data types")
count(DE_data, data_type)

# data_type cluster contrast seqnames     start       end  width strand source type  gene_id         gene_version gene_name       gene_type      vlmf_logFC vlmf_AveExpr vlmf_t vlmf_P.Value vlmf_adj.P.Val vlmf_B cell_type_broad
# <chr>     <fct>   <chr>    <fct>        <int>     <int>  <int> <fct>  <fct>  <fct> <chr>           <chr>        <chr>           <chr>               <dbl>        <dbl>  <dbl>        <dbl>          <dbl>  <dbl> <fct>          
# 1 sn_fine   Astro.1 carrier  chr15     82339993  82349475   9483 -      HAVANA gene  ENSG00000278662 6            GOLGA6L10       protein_coding     -1.55          1.83  -6.46  0.000000659         0.0111   5.60 Astro          
# 2 sn_fine   Astro.1 carrier  chr22     38818452  38844028  25577 -      HAVANA gene  ENSG00000221890 6            NPTXR           protein_coding      0.751         4.29   5.73  0.00000436          0.0310   4.28 Astro          
# 3 sn_fine   Astro.1 carrier  chr1      99264292  99309590  45299 +      HAVANA gene  ENSG00000117600 14           PLPPR4          protein_coding      1.50          2.27   5.64  0.00000552          0.0310   3.76 Astro          
# 4 sn_fine   Astro.1 carrier  chr6      40344346  40387194  42849 -      HAVANA gene  ENSG00000226070 2            LINC00951       lncRNA              1.19          1.81   5.47  0.00000885          0.0373   3.21 Astro  


key_genes <- read_csv(here(data_dir, "key_genes_context.csv"))
# gene    cell_type mediator_role data_type      cluster_key evidence_tier note                                                
# <chr>   <chr>     <chr>         <chr>          <chr>       <chr>         <chr>                                               
#     1 NPTXR   Astro.1   source        cell_type_anno Astro.1     top           mediator source; validates pooled + L5/L2.3         
# 2 NPTXR   Oligo.3   target        cell_type_anno Oligo.3     top           mediator target; validates pooled + WM/L5/L1b/Vasc  
# 3 FZD8    Astro.2   source        SpX            Astro.2_L6  top           strongest point on WM-proximal gradient             
# 4 CABLES1 Oligo.3   highlight     Oligo_3_Astro  nn_Astro.3  top           masked when pooled; only emerges near astrocytes    
# 5 CPNE4   Oligo.3   highlight     SpX            Oligo.3_L6  top           Match_all; also independently sig. in Visium L6     
# 6 CD9     Oligo.3   highlight     Visium         L6~Sp09D04  second        myelin gene; off Xenium panel; Visium-only evidence 
# 7 NRG1    Oligo.5   hypothesis    MOFA           NA          hypothesis    top MOFA loading; ERBB3 ligand; not yet DE-validated


#### Analysis context ####

## process one cell type / spatial domain at a time - change this and rerun
## to build the panel for a different context
analysis_context <- "Oligo.3"

key_genes_context <- key_genes |>
    filter(cell_type == analysis_context)

gene_levels <- key_genes_context$gene

#### Filter to key genes, discovery-gated Xenium validation ####

## a single column to identify "which cell type does this row's test belong
## to" regardless of data_type - Xenium_cell_type_anno has no cell_type_anno
## column of its own (cluster IS the cell type there), while Xenium_SpX and
## Xenium_Oligo.3_Astro carry a genuine cell_type_anno distinct from their
## (compound, stratified) cluster names
DE_data <- DE_data |>
    mutate(cell_type_match = coalesce(cell_type_anno, as.character(cluster)))

## sn_fine hits specifically in this context define "discovered" - strictly
## restricted to the analysis_context cluster, not any other sn_fine cluster
discovery_hits <- DE_data |>
    filter(data_type == "sn_fine", cluster == analysis_context,
           gene_name %in% gene_levels, vlmf_adj.P.Val < 0.05) |>
    distinct(gene_name)

## discovery rows: sn_fine filtered down to just this context's cluster;
## Visium kept at whatever SpD it was significant in (SpDs aren't named the
## same as cell types, so there's no cluster name to match on there)
discovery_key <- DE_data |> 
    filter(
        (data_type == "sn_fine" & cluster == analysis_context) | data_type == "Visium",
        gene_name %in% gene_levels, vlmf_adj.P.Val < 0.05
    ) |>
    mutate(context = paste0(data_type, "|", cluster), context_type = "discovery")

## Xenium_cell_type_anno is already one row per gene x cell type ("pooled") -
## restrict to this context's cell type + discovered genes + nominal cutoff
xenium_pooled <- DE_data |>
    filter(data_type == "Xenium_cell_type_anno", cell_type_match == analysis_context,
           gene_name %in% discovery_hits$gene_name, vlmf_P.Value < 0.1) |>
    mutate(context = paste0(data_type, "|", cluster), context_type = "xenium_pooled")

## Xenium_Oligo.3_Astro and Xenium_SpX are stratified into many sub-contexts
## within this cell type (astrocyte proximity, spatial domain) - restrict to
## this context's cell type + discovered genes, then collapse to the single
## best (lowest p) context per gene per data type
best_stratified_context <- function(dt) {
    DE_data |>
        filter(data_type == dt, cell_type_match == analysis_context,
               gene_name %in% discovery_hits$gene_name, vlmf_P.Value < 0.1) |>
        group_by(gene_name) |>
        # slice_min(vlmf_P.Value, n = 1, with_ties = FALSE) |>
        ungroup() |>
        mutate(context = paste0(data_type, "|", cluster), context_type = "xenium_best")
}

xenium_best <- c("Xenium_Oligo.3_Astro", "Xenium_SpX") |>
    map(best_stratified_context) |>
    list_rbind()

DE_data_key <- bind_rows(discovery_key, xenium_pooled, xenium_best) |>
    mutate(gene_name = factor(gene_name, levels = gene_levels))

message(n_distinct(DE_data_key$context), " contexts across ", n_distinct(DE_data_key$gene_name), " genes in the ", analysis_context, " context")
count(DE_data_key, gene_name, data_type, context_type)

#### Plot ####

p_de <- DE_data_key |> 
    ggplot(aes(x = context, y = gene_name, fill = vlmf_t)) +
    geom_tile() +
    facet_wrap(~data_type, nrow = 1, scales = "free_x") +
    scale_fill_gradient2(low = APOE_carrier_colors[["E2+"]], 
                         high = APOE_carrier_colors[["E4+"]], 
                         midpoint = 0) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
    

p_de

ggsave(p_de, filename = here(plot_dir, sprintf("%s_contexts_heatmap.pdf", analysis_context)), width = 8, height = 4)


#### Session info ####

session_info()

