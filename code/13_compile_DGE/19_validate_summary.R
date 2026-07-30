## Louise Huuki-Myers July 2026
## Check total sets of 

#### set up ####
library("here")
library("tidyverse")
library("qs2")
library("sessioninfo")
library("patchwork")
library("ComplexHeatmap")
library("circlize")

data_dir <- here("processed-data", "13_compile_DGE", "19_validate_summary")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "19_validate_summary")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
# APOE_carrier_colors
# E2+       E4+
#"#51B8B1" "#D97D59"

## CHECK: logFC_Heatmap() looks up `cluster_levels` from the global env for
## row ordering (it's not a function argument) - loaded per-datatype below,
## right before each call, so it's always the current resolution's levels.

## fill_stat-adaptable logFC/t-stat heatmap helper (logFC_Heatmap / logFC_Heatmap_contrast)
source(here("code", "13_compile_DGE", "logFC_heatmap.R")) 


#### Validation data ####

## Xenium spatial validation layers, brought in as the 4th resolution.
## `_wSN` files already carry a `validate` column (direction-concordant +
## meets the nominal validation threshold vs. discovery) - use that directly
## rather than re-deriving significance from a raw p-value/FDR column.
validation_data_types <- c("Xenium_cell_type_anno", "Xenium_cell_type_anno_SpX", "Xenium_Oligo.3_Astro")

load_DE_valid_data <- function(datatype){
    DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s_wSN.Rds", datatype))
    stopifnot("DE validation data file not found" = file.exists(DE_data_fn))
    message("loading ", datatype, ": ", DE_data_fn)
    readRDS(DE_data_fn)
}

DE_valid_data <- validation_data_types |>
    map(load_DE_valid_data) |>
    list_rbind() |>
    mutate(cluster = ifelse(data_type == "Xenium_cell_type_anno_SpX", cluster_SpX, as.character(cluster)))

all(anchor_genes %in% DE_valid_data$gene_name)

message(nrow(DE_valid_data), " rows across ", n_distinct(DE_valid_data$data_type), " validation data types")
count(DE_valid_data, data_type)

# quick look at which clusters actually validate, per resolution
DE_valid_data |>
    filter(validate) |> 
    count(data_type, cluster)
    
valid_genes_summary <- DE_valid_data |>
    filter(validate) |>
    group_by(data_type, gene_name, cell_type_anno) |>
    summarise(n = n(), enviro = paste0(unique(cluster), collapse = ", ")) |>
    arrange(-n)  |>
    mutate(enviro2=ifelse(n>1, "multi", enviro))

valid_genes_summary |> 
    filter(cell_type_anno == "Oligo.3") |> 
    arrange(gene_name) |> 
    print(n = 53) 

DE_valid_data |> filter(data_type == "Xenium_SpX") |> distinct(cluster)


#### bar plot ####


DE_valid_data |> 
    filter(cell_type_anno == "Oligo.3", validate) |>
    ggplot(aes(x = data_type)) +
    geom_bar(aes(fill = cluster))

valid_genes_summary |>
    filter(cell_type_anno == "Oligo.3") |>
    mutate(enviro2 = ifelse(n > 1, "multi", enviro)) |>
    ggplot(aes(x = data_type)) +
    geom_bar(aes(fill = enviro2))

# valid_genes_summary |>
#     filter(cell_type_anno=="Oligo.3") |> 
#     write_csv(here(data_dir, "valid_genes_summary_test.csv"))
    
valid_genes_summary_bar_text <- valid_genes_summary |>
    filter(cell_type_anno=="Oligo.3") |>
    mutate(enviro2=fct_relevel(factor(enviro2), "multi")) |>
    arrange(data_type, enviro2) |>
    group_by(data_type) |>
    mutate(y_pos=row_number() - 0.5) |>
    ungroup() |>
    ggplot(aes(x=data_type)) +
    geom_bar(aes(fill=enviro2)) +
    geom_text(aes(y=y_pos, label=gene_name), size=3, hjust=0.5) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1))

dt = "Xenium_cell_type_anno"
lookup <- load_cluster_lookup(dt)
cluster_colors <- lookup$cluster_colors  # logFC_Heatmap() auto-builds a row/col color annotation from this
cluster_levels <- lookup$cluster_levels  # logFC_Heatmap() reads this from the global env

dt_data <- DE_valid_data |> filter(data_type == dt)

## match the clusterd order from discovery dataset
anchor_genes2 <- c("NPTXR", "CPNE4", "CABLES1", "FZD8", "FOS", "STAT4", "TLR2", "STAT1", "SOX10", "MBP", "MAL", "MAG", "KLK6", "OPALIN", "PLP1")

## full version - every cluster in this resolution
logFC_Heatmap(
    data = dt_data,
    gene_list = anchor_genes2,
    title = "anchor_genes",
    datatype = dt,
    fill_stat = "vlmf_xenium_t",
    flip = TRUE,
    order_genes = FALSE,  # keep the fixed thematic gene order above, not a data-driven one
    save = TRUE,
    h = 5, 
    w = (length(unique(dt_data$cluster))*0.3) + 2,
    valid_anno = "validate"
)

## filtered version - only clusters validation
valid_clusters <- dt_data |>
    filter(gene_name %in% anchor_genes2, validate) |>
    distinct(cluster) |>
    pull(cluster)

if(length(valid_clusters) == 0){
    message(dt, ": no validation clusters - skipping filtered heatmap")
} else {
    logFC_Heatmap(
        data = dt_data |> filter(cluster %in% valid_clusters), ## use sig_clusters to include Inhib.Vip
        gene_list = anchor_genes2,
        title = "anchor_genes_sig_clusters",
        datatype = dt,
        fill_stat = "vlmf_xenium_t",
        flip = TRUE,
        order_genes = FALSE,
        save = TRUE,
        h = 5, 
        w = length(sig_clusters) *0.9 ,
        valid_anno = "validate"
    )
}

#### Xenium SpX validation heatmap ####
## Will have to provide custom annotations for Xenium SpX
dt = "Xenium_cell_type_anno_SpX"

# lookup <- load_cluster_lookup(dt)
# cluster_colors <- lookup$cluster_colors  # logFC_Heatmap() auto-builds a row/col color annotation from this
# cluster_levels <- lookup$cluster_levels  # logFC_Heatmap() reads this from the global env

dt_data <- DE_valid_data |> 
    filter(data_type == dt) |> 
    mutate(cluster = cluster_SpX) 

dt_data |> 
    filter(gene_name == "CABLES1", cell_type_anno == "Oligo.3") |>
    select(gene_name, cluster, cell_type_anno, SpX, vlmf_xenium_P.Value, vlmf_xenium_t, validate)

# gene_name cluster       cell_type_anno SpX        vlmf_xenium_P.Value vlmf_xenium_t validate
# 6 CABLES1   Oligo.3_L6    Oligo.3        L6~SpX9                 0.228          1.25  FALSE   ????
# 7 CABLES1   Oligo.3_Vasc  Oligo.3        Vasc~SpX3               0.340         -0.977 FALSE   
# 8 CABLES1   Oligo.3_WM    Oligo.3        WM~SpX2                 0.0512         2.07  TRUE    
# 9 CABLES1   Oligo.3_WMtz  Oligo.3        WMtz~SpX8               0.0931         1.77  TRUE

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpX_colors.Rdata"), verbose = TRUE)

dt_xenium_spx_anno <- dt_data |>
    distinct(cluster, cell_type_anno, SpX) |>
    column_to_rownames("cluster")

ha_xenium_spx <- HeatmapAnnotation(df = dt_xenium_spx_anno,
                                   col = list(cell_type_anno = cell_type_colors$anno,
                                              SpX = SpX_colors))

# mutate(SpX_simple =str_remove(SpX, "~.*"))

## full version - every cluster in this resolution
logFC_Heatmap(
    data = dt_data,
    gene_list = anchor_genes2,
    title = "anchor_genes",
    datatype = dt,
    fill_stat = "vlmf_xenium_t",
    flip = TRUE,
    col_anno = ha_xenium_spx, 
    order_genes = FALSE,  # keep the fixed thematic gene order above, not a data-driven one
    save = TRUE,
    h = 5, 
    w = (length(unique(dt_data$cluster))*0.3) + 2,
    valid_anno = "validate"
)

## filtered version - only clusters validation
valid_clusters <- dt_data |>
    filter(gene_name %in% anchor_genes2, validate) |>
    distinct(cluster) |>
    pull(cluster)

if(length(valid_clusters) == 0){
    message(dt, ": no validation clusters - skipping filtered heatmap")
} else {
    
    dt_xenium_spx_anno <- dt_data |> 
        filter(cluster %in% valid_clusters) |>
        distinct(cluster, cell_type_anno, SpX) |>
        column_to_rownames("cluster")
    
    ha_xenium_spx <- HeatmapAnnotation(df = dt_xenium_spx_anno,
                                       col = list(cell_type_anno = cell_type_colors$anno[unique(dt_xenium_spx_anno$cell_type_anno)],
                                                  SpX = SpX_colors[unique(dt_xenium_spx_anno$SpX)]))
    
    logFC_Heatmap(
        data = dt_data |> filter(cluster %in% valid_clusters),
        gene_list = anchor_genes2,
        title = "anchor_genes_sig_clusters",
        datatype = dt,
        col_anno = ha_xenium_spx, 
        fill_stat = "vlmf_xenium_t",
        flip = TRUE,
        order_genes = FALSE,
        save = TRUE,
        h = 5, 
        w = length(valid_clusters) *0.9 ,
        valid_anno = "validate"
    )
}




#### Session info ####

session_info()
