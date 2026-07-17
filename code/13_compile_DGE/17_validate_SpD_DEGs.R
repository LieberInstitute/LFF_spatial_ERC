## Louise Huuki-Myers July 2026
## Load carrier DGE results across data types (Visium SpDs + Xenium SpX/cell type angles) and combine for cross-platform validation.

#### set up ####
library("here")
library("tidyverse")
library("qs2")
library("sessioninfo")
library("ggrepel")
library("ComplexHeatmap")


data_dir <- here("processed-data", "13_compile_DGE", "17_validate_SpD_DEGs")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "17_validate_SpD_DEGs")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))

#### load data ### 

DE_data_visium <- readRDS(here("processed-data","13_compile_DGE","01_compile_DGE","Visium" ,"DGE_results_carrier_Visium.Rds")) |>
    mutate(SpD_simple =str_remove(cluster, "~.*"))

DE_data_xenium <- readRDS(here("processed-data","13_compile_DGE","01_compile_DGE","Xenium_SpX" ,"DGE_results_carrier_Xenium_SpX_wSN.Rds")) |>
    mutate(SpX_simple =str_remove(SpX, "~.*"))

DE_data_visium |> distinct(SpD_simple)
DE_data_xenium |> distinct(SpX_simple)

Sp_lookup <- tribble(
    ~SpX_simple, ~SpD_simple,
    "Inhib",     "Inhib",
    "L1a",       "L1",
    "L1b",       "L1",
    "L2.3",      "L2.3",
    "L5",        "L5",
    "L6",        "L6",
    "Vasc",      "Vasc",
    "WM",        "WM",
    "WMtz",      "WM.uf"
)
    

DE_data_visium |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster)
DE_data_visium |> filter(vlmf_adj.P.Val < 0.05, gene_name %in% DE_data_xenium$gene_name) |> count(cluster)


DE_data_combine <- DE_data_visium |> 
    select(cluster, SpD_simple, gene_name, starts_with("vlmf")) |>
    rename_with(~ str_replace(.x, "vlmf_", "vlmf_visium_"), starts_with("vlmf_")) |>
    left_join(Sp_lookup, by = join_by(SpD_simple), relationship = "many-to-many") |>
    inner_join(DE_data_xenium |> rename(Xenium_cluster = cluster),
              relationship = "many-to-many") |>
    mutate(signif_visium = vlmf_visium_adj.P.Val < 0.05,
           dir_match_vis = (vlmf_xenium_t > 0) == (vlmf_visium_t > 0),
           match_class = case_when(signif_xenium & dir_match & signif_sn & dir_match_vis & signif_visium~ "Match_all",
                                   signif_xenium & signif_visium & dir_match_vis ~ "Match_Visium",
                                   signif_xenium & dir_match & signif_sn ~ "Match_sn",
                                   TRUE ~ "No Match"))


DE_data_combine |> 
    count(match_class)

DE_data_combine |> 
    filter(visium_sig, signif_xenium, dir_match_vis) |>
    count(gene_name)

# gene_name     n
# <chr>     <int>
# 1 CNDP1         4
# 2 CNTN2         2
# 3 CPNE4        13
# 4 KLK6          2
# 5 MAG          12
# 6 MAL           2
# 7 MOG           6
# 8 MYRF          5
# 9 NTNG2        12
# 10 NWD2          4
# 11 UGT8          2


DE_data_combine |> 
    filter(match_class == "Match_all") |>
    select(gene_name, SpD_simple, SpX_simple, cell_type_anno, vlmf_visium_t, vlmf_visium_adj.P.Val , vlmf_sn_t, vlmf_sn_adj.P.Val , vlmf_xenium_t, vlmf_xenium_adj.P.Val, dir_match, dir_match_vis)
# gene_name SpD_simple SpX_simple cell_type_anno vlmf_visium_t vlmf_visium_adj.P.Val vlmf_sn_t vlmf_sn_adj.P.Val vlmf_xenium_t vlmf_xenium_adj.P.Val dir_match dir_match_vis
# <chr>     <chr>      <chr>      <fct>                  <dbl>                 <dbl>     <dbl>             <dbl>         <dbl>                 <dbl> <lgl>     <lgl>        
#1 CPNE4     L6         L6         Oligo.3                 5.06                0.0429      3.39            0.0350          2.02                 0.543 TRUE      TRUE 

DE_data_combine |>
    filter(gene_name == "CPNE4")|>
    select(gene_name, SpD_simple, SpX_simple, cell_type_anno, vlmf_visium_t, vlmf_visium_adj.P.Val , vlmf_sn_t, vlmf_sn_adj.P.Val , vlmf_xenium_t, vlmf_xenium_P.Value, dir_match, dir_match_vis)

visium_t <- DE_data_visium |>
    filter(gene_name == "KLK6")|>
    select(gene_name, SpD_simple, vlmf_t, vlmf_adj.P.Val) |>
    filter(vlmf_adj.P.Val < 0.05) |>
    pull(vlmf_t)


DE_data_combine_hits <- DE_data_combine |>
    filter(signif_visium & signif_xenium & dir_match_vis) |>
    select(gene_name, SpD_simple, SpX_simple, cell_type_anno, 
           vlmf_visium_t, vlmf_visium_adj.P.Val , vlmf_sn_t, vlmf_sn_adj.P.Val,
           vlmf_sn_AveExpr , vlmf_xenium_t, 
           vlmf_xenium_P.Value, vlmf_xenium_AveExpr, 
           dir_match, dir_match_vis, match_class)


DE_data_combine_hits |> group_by(SpD_simple, gene_name) |> summarize(visium = length(unique(SpD_simple)), xenium = sum(vlmf_xenium_P.Value < 0.1 & dir_match_vis))


valid_genes<- DE_data_combine_hits |>
    distinct(gene_name) |>
    pull(gene_name)

valid_SpD <- DE_data_combine_hits |>
    distinct(SpD_simple) |>
    pull(SpD_simple)


valid_pairs <- DE_data_combine_hits |>
    distinct(gene_name, SpX_simple, vlmf_visium_t, vlmf_visium_adj.P.Val) 

match_class_colors <- c("Match_Visium" = "#a183ff",
                        "Match_all" = "#d20057",
                        "Match_sn" = "#00a293",
                        "No Match" = "grey")

SpD_valid_scatter <-map(valid_SpD, ~DE_data_combine |>
                            filter(signif_visium, SpD_simple == .x)|>
                            # filter(gene_name %in% valid_genes)|>
                            ggplot(aes(x = vlmf_sn_t, y = vlmf_xenium_t, color = match_class)) +
                            geom_hline(yintercept = 0, color = "grey30")+
                            geom_vline(xintercept = 0, color = "grey30") +
                            geom_point() +
                            geom_text_repel(aes(label = cell_type_anno)) +
                            geom_hline(aes(yintercept = vlmf_visium_t), linetype = "dashed", color = "red") +
                            geom_vline(aes(xintercept = vlmf_visium_t), linetype = "dashed", color = "red") +
                            facet_grid(SpD_simple~gene_name) +
                            scale_color_manual(values = match_class_colors) +
                            theme_bw() 
)

walk2(SpD_valid_scatter, valid_SpD, ~ggsave(.x, filename = here(plot_dir, sprintf("SpD_valid_t_scatter_%s.png", .y)), width = 9, height = 5))

# DE_data_combine |>
#     filter(signif_visium, SpD_simple == "L6")|>
#     # filter(gene_name %in% valid_genes)|>
#     ggplot(aes(x = vlmf_sn_t, y = vlmf_xenium_t, color = match_class)) +
#     geom_point() +
#     geom_text_repel(aes(label = cell_type_anno)) +
#     geom_hline(aes(yintercept = vlmf_visium_t), linetype = "dashed", color = "red") +
#     facet_grid(SpD_simple~gene_name) +
#     theme_bw() +
#     geom_hline(yintercept = 0, color = "grey30")+
#     geom_vline(xintercept = 0, color = "grey30")
#     

DE_data_combine_hits |> filter(gene_name == "CPNE4", vlmf_xenium_P.Value < 0.1) |> arrange(-vlmf_xenium_t)
DE_data_combine_hits |> filter(gene_name == "UGT8", vlmf_xenium_P.Value < 0.1) |> arrange(-vlmf_xenium_t)
DE_data_combine_hits |> filter(gene_name == "KLK6", vlmf_xenium_P.Value < 0.1) |> arrange(-vlmf_xenium_t)


DE_data_combine_hits |>
    write_csv(here(data_dir, "SpD_DEG_xenium_valid.csv"))


DE_data_xenium |>  select(gene_name,SpX, cell_type_anno, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_xenium_AveExpr)

    
    
DE_data_xenium |> 
    filter(gene_name == "CPNE4") |>
    ggplot(aes(x = cell_type_anno, y = SpX, fill = vlmf_xenium_t)) +
    geom_tile()

#### Volcano plots ####
valid_pairs 

map(unique(valid_pairs$SpX_simple), function(s){
    
    pairs_s <- valid_pairs |> 
        filter(SpX_simple == s) |> 
        mutate(t_anno = sprintf("visium t = %s", round(vlmf_visium_t, 2)))
    
    visium_validation_volcano <- DE_data_combine |> 
        filter(gene_name %in% pairs_s$gene_name, SpX_simple == s) |>
        ggplot(aes(x = vlmf_xenium_logFC, y = -log10(vlmf_xenium_P.Value), color = match_class)) +
        geom_point() +
        geom_text_repel(aes(label = ifelse(match_class != "No Match", as.character(cell_type_anno), ""))) +
        scale_color_manual(values = match_class_colors) +
        geom_label(
            data = pairs_s,
            aes(x = Inf, y = -Inf, label = t_anno),
            size = 4, vjust = "inward", hjust = "inward",
            color = "black"
        ) +
        facet_grid(SpX_simple~gene_name) +
        theme_bw() +
        theme(legend.position = "bottom")
    
    ggsave(visium_validation_volcano, filename = here(plot_dir, sprintf("visium_validation_volcano_%s.png", s)), width = nrow(pairs_s)*2 + 2, height = 5)
    
})

DE_data_combine |> 
    filter(gene_name == "CPNE4", SpX_simple == "L6") |>
    ggplot(aes(x = vlmf_xenium_logFC, y = -log10(vlmf_xenium_P.Value), color = match_class)) +
    geom_point() +
    geom_text_repel(aes(label = cell_type_anno)) +
    scale_color_manual(values = match_class_colors) +
    facet_grid(SpX_simple~gene_name) +
    theme_bw()


DE_data_combine |> 
    filter(gene_name %in% valid_genes) |>
    ggplot(aes(x = vlmf_xenium_logFC, y = -log10(vlmf_xenium_P.Value), color = match_class)) +
    geom_point() +
    geom_text_repel(aes(label = ifelse(match_class != "No Match", as.character(cell_type_anno), ""))) +
    scale_color_manual(values = match_class_colors) +
    facet_grid(SpX_simple~gene_name) +
    theme_bw()

sig_facets <- DE_data_combine |>
    filter(gene_name %in% valid_genes) |>
    group_by(SpX_simple, gene_name) |>
    summarise(is_sig = any(vlmf_visium_adj.P.Val < 0.05, na.rm = TRUE),
              # de = case_when(is_sig & vlmf_visium_t ~ "UP",
              #                is_sig & vlmf_visium_t ~ "DOWN",
              #                TRUE ~"None"),
              .groups = "drop")

sig_facets |> filter(is_sig)

DE_data_combine |>
    filter(gene_name %in% valid_genes) |>
    ggplot(aes(x = vlmf_xenium_logFC, y = -log10(vlmf_xenium_P.Value), color = match_class)) +
    geom_rect(data = sig_facets, aes(fill = is_sig),
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
              inherit.aes = FALSE, alpha = 0.15) +
    scale_fill_manual(values = c(`TRUE` = "gold", `FALSE` = NA), guide = "none") +
    geom_point() +
    geom_text_repel(aes(label = ifelse(match_class != "No Match", as.character(cell_type_anno), ""))) +
    scale_color_manual(values = match_class_colors) +
    facet_grid(SpX_simple~gene_name) +
    theme_bw()


#### Heatmaps ####

source("SpX_tstat_heatmap.R")
source("SpX_tstat_heatmap2.R")



pdf(here(plot_dir, "Visium_validation_SpX_heatmap.pdf"))
walk(valid_genes, ~SpX_tstat_heatmap(DE_data_xenium, gene=.x, datatype="Xenium_SpX", save=FALSE, order_celltypes = FALSE, visium_data = DE_data_visium))
dev.off()



## TODO debug - actually its okay - just missing data in snRNA-seq 
DE_data_visium |> filter(gene_name == "NWD2") |> arrange(vlmf_P.Value) 
DE_data_xenium |> filter(gene_name == "NWD2", signif_xenium, grepl("L6", SpX)) |> arrange(vlmf_xenium_adj.P.Val) 

DE_data_combine_hits |> filter(gene_name == "NWD2") |> arrange(vlmf_xenium_P.Value) 
DE_data_combine_hits |> filter(gene_name == "CAPN3") |> arrange(vlmf_xenium_P.Value) 
    