## Louise Huuki-Myers, May 2026
## Compile and plot all DGE data from Xenium

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("GGally")
library("getopt")
library("broom")
library("ComplexHeatmap")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt$datatype  <- "Xenium_cell_type_anno_SpX"

data_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### colors and factors ####
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpX_colors.Rdata"), verbose = TRUE)
SpX_levels <- names(SpX_colors)

cell_type_broad_levels <- names(cell_type_colors$broad)
cell_type_broad_levels <- cell_type_broad_levels[cell_type_broad_levels != "Other"]

if(opt$datatype == "Xenium_cell_type_anno_SpX"){
    
    cluster_levels <- names(cell_type_colors$anno)
    
    load(here("processed-data", "SpX_colors.Rdata"), verbose = TRUE)
    SpX_levels <- names(SpX_colors)
    
    SpX_levels_simple <- gsub("~.*", "", SpX_levels)
    
    SpX_colors_simple <- SpX_colors
    names(SpX_colors_simple) <- SpX_levels_simple

        
    cell_type_anno_SpX_level_tb <- expand_grid(cell_type_anno = cell_type_levels, SpX_simple = SpX_levels_simple) |>
        mutate(cell_type_anno_SpX = paste0(cell_type_anno, '_', SpX_simple))
    
    cluster_spx_levels <- cell_type_anno_SpX_level_tb$cell_type_anno_SpX
    
} else if(opt$datatype == "Xenium_Oligo.3_Astro_SpX"){
    
    cluster_levels <- c("nnA_far_APOE_low", "nnA_far_APOE_high", "nnA_near_APOE_low", "nnA_near_APOE_high")
    
    cluster_colors <- c("nnA_far_APOE_low" = "#87219a",
                        "nnA_far_APOE_high" = "#487800",
                        "nnA_near_APOE_low" = "#b0b3ff",
                        "nnA_near_APOE_high" = "#01d9b4"
                        # "doublet_Oligo3_Astro" = "#ef2396"
    )
}


#### voomLmFit data ####
vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "06_Clusterwise_voomLmFit_Xenium", paste0("vlmf_", opt$datatype)),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub(sprintf("voomLmFit_%s_|.rds", opt$datatype), "", basename(.x)))

## read data
vlmf_data <- map(vlmf_fn, readRDS)

# vlmf_data <- list_transpose(vlmf_data)

names(vlmf_data)
head(vlmf_data[[1]])

# all(names(vlmf_data) %in% cluster_spx_levels)
# cluster_spx_levels <- cluster_levels[cluster_spx_levels %in% names(vlmf_data)]

vlmf_data_tb <- do.call("rbind", vlmf_data) |>
    dplyr::rename(vlmf_logFC = logFC,
                  vlmf_AveExpr = AveExpr,
                  vlmf_t = t,
                  vlmf_P.Value = P.Value,
                  vlmf_adj.P.Val = adj.P.Val,
                  vlmf_B = B )  |>
    mutate(cluster_SpX = cluster,
           cluster = droplevels(factor(gsub("_[^_]+$", "", cluster), cluster_levels)),
           SpX = factor(SpX, names(SpX_colors)),
           cell_type_anno = droplevels(factor(cell_type_anno, names(cell_type_colors$anno)))) |>
    as_tibble() 

cluster_levels <- levels(vlmf_data_tb$cluster)
levels(vlmf_data_tb$SpX)

vlmf_data_tb |> filter(vlmf_P.Value < 0.05) |> count(SpX) 
vlmf_data_tb |> filter(vlmf_P.Value < 0.05) |> count(cluster_SpX) |> arrange(-n) |> print(n = 35)

vlmf_data_tb |> count(SpX)

vlmf_model_summary <- vlmf_data_tb |> 
                                   mutate(mod = "carrier") |>
                                   group_by(SpX, cluster, mod) |>
                                   summarize(n_genes= n(),
                                             n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
                                             n_pval10 = sum(vlmf_P.Value < 0.05),
                                             nUP = sum(vlmf_P.Value < 0.05 & vlmf_logFC > 0),
                                             nDown = sum(vlmf_P.Value < 0.05 & vlmf_logFC < 0))

vlmf_model_summary |> filter(n_FDR05 > 0)
vlmf_model_summary |> arrange(-n_pval10)
# vlmf_model_summary |> filter(cluster == "Oligo.3")

## n signif bar plots
vlmf_model_summary_bar <- vlmf_model_summary |>
    ggplot(aes(x = SpX, y = n_pval10, fill = SpX)) +
    geom_col() +
    geom_text(aes(label = n_pval10), vjust=-.5) +
    scale_fill_manual(values = SpX_colors) +
    facet_wrap(~cluster) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "p-value < 0.10") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None")

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar.png", opt$datatype)), height = 12, width = 12)

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_pval10, -n_genes, -n_FDR05) |>
    mutate(nDown = -1*nDown) |>
    pivot_longer(!c(cluster, SpX, mod), names_to = "reg", values_to = "n_genes") |>
    # filter(startsWith(mod, "apoe") | mod %in% c("E4E4", "carrier")) |>
    filter(mod == "carrier") |>
    ggplot(aes(x = SpX, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = abs(n_genes))) +
    facet_wrap(~cluster) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "p-value < 0.10") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar_reg.png", opt$datatype)), width = 10)

## save summary
write.csv(vlmf_model_summary, file = here(data_dir, sprintf("vlmf_model_summary_%s.csv", opt$datatype)))

#### vlmf volcano plots ####

custom_volcano <- function(data, Pval_cut = 0.10, model_name){
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("both", paste("Pval<", Pval_cut) , "abs(logFC)>1" )
    
    volcano <- data |>
        mutate(DE_class = case_when(vlmf_P.Value < Pval_cut & abs(vlmf_logFC) > 1  ~ "both",
                         vlmf_P.Value < Pval_cut ~ paste("Pval<", Pval_cut),
                         abs(vlmf_logFC) > 1 ~ "abs(logFC)>1",
                         TRUE ~ "None")) |>
        ggplot(aes(x = vlmf_logFC, y = -log10(vlmf_P.Value), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        scale_color_manual(values = signif_colors) +
        # facet_wrap(~cluster) +
        facet_wrap(~SpX) +
        theme_bw() +
        labs(title = model_name)
    
    if(nrow(data) > 1000){
        volcano <- volcano + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) 
    } else {
        volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
    }
    ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_plot_%s.png", model_name)), height = 10, width = 10)
}


custom_volcano_ct_fine <- function(data, mod_name = "carrier"){
    
    map(cluster_levels, 
        ~data |> 
            filter(cluster == .x) |>
            custom_volcano(model_name = paste0(mod_name,"-", .x)) )
    
}


## plot volcanos
custom_volcano_ct_fine(vlmf_data_tb)


#### Save vlmf data & Add sn results ####
sn_DEG_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds")) |>
    dplyr::select(cell_type_anno = cluster, gene_id, gene_name, starts_with("vlmf")) |>
    rename_with(~ str_replace(.x, "vlmf_", "vlmf_sn_"), starts_with("vlmf_")) 

sn_DEG_data |> count(vlmf_sn_adj.P.Val < 0.05, gene_name %in% vlmf_data_tb$gene_name)
# sn_DEG_data |> count(gene_name %in% vlmf_data_tb$gene_name)

if(opt$datatype == "Xenium_cell_type_anno_SpX"){
    
    vlmf_data_tb_xenium <- vlmf_data_tb |>
        rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_"))|> 
        inner_join(sn_DEG_data) |>
        mutate(signif_xenium = vlmf_xenium_P.Value < 0.1, 
               signif_sn = vlmf_sn_adj.P.Val < 0.05,
               signif_both = signif_xenium & signif_sn,
               dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0),
               validate = signif_both & dir_match)
    
    
} else if(opt$datatype == "Xenium_Oligo.3_Astro_SpX"){
    
    sn_DEG_data <- sn_DEG_data |>
        filter(cell_type_anno == "Oligo.3") 
    
    vlmf_data_tb_xenium <- vlmf_data_tb |>
        rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_"))|> 
        inner_join(sn_DEG_data) |>
        mutate(signif_xenium = vlmf_xenium_P.Value < 0.1, 
               signif_sn = vlmf_sn_adj.P.Val < 0.05,
               signif_both = signif_xenium & signif_sn,
               dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0),
               validate = signif_both & dir_match)
}



vlmf_data_tb_xenium |> count(cluster_SpX)

## save data
saveRDS(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_carrier_%s.Rds", opt$datatype)))
write.csv(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_carrier_%s.csv", opt$datatype)), row.names = FALSE)

## save data
saveRDS(vlmf_data_tb_xenium, file = here(data_dir, sprintf("DGE_results_carrier_%s_wSN.Rds", opt$datatype)))
write.csv(vlmf_data_tb_xenium, file = here(data_dir, sprintf("DGE_results_carrier_%s_wSN.csv", opt$datatype)), row.names = FALSE)

#### Summarize validation experiments ####

validation_summary <- vlmf_data_tb_xenium |> 
    filter(!is.na(vlmf_sn_P.Value)) |> 
    group_by(cluster, SpX) |>
    summarise(n_signif_xenium = sum(signif_xenium),
              n_signif_sn = sum(signif_sn),
              n_signif_both = sum(signif_both),
              n_validate = sum(signif_both & dir_match),
              percent_valid = 100*sum(validate)/n_signif_sn)

validation_summary |>
    filter(n_signif_sn > 0)

validation_summary |> arrange(-n_validate)

validation_summary |>
    filter(cluster == "Oligo.3")

vlmf_data_tb_xenium |>
    filter(validate) |> 
    select(cluster, SpX, gene_name, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_sn_t) |> 
    arrange(cluster, SpX) |>
    print(n = 79)

vlmf_data_tb_xenium |> 
    filter(validate) |> 
    group_by(gene_name) |>
    count() |> 
    arrange(-n) 


vlmf_data_tb_xenium |> 
    # filter(gene_name == "IGFBP5") |>
    filter(gene_name == "CABLES1", validate) |>
    # filter(SpX == "L1b~SpX7", signif_both) |>
    # filter(SpX == "L6~SpX9", signif_both) |>
    # filter(SpX == "WMtz~SpX8", signif_both) |>
    select(cluster_SpX, gene_name, vlmf_sn_t, vlmf_sn_adj.P.Val, vlmf_xenium_t, vlmf_xenium_logFC, vlmf_xenium_P.Value, dir_match,validate) |>
    arrange(-vlmf_xenium_logFC)

# cluster gene_id     gene_name vlmf_sn_logFC vlmf_sn_adj.P.Val vlmf_xenium_logFC vlmf_xenium_P.Value dir_match validate
# <fct>   <chr>       <chr>             <dbl>             <dbl>             <dbl>               <dbl> <lgl>     <lgl>
# 1 Oligo.3 ENSG000001… SLC17A7           1.47             0.0388             0.234              0.0453 TRUE      TRUE
# 2 Oligo.3 ENSG000001… ENC1              1.61             0.0238             0.246              0.0530 TRUE      TRUE
# 3 Oligo.3 ENSG000002… NPTXR             1.21             0.0178             0.399              0.0497 TRUE      TRUE
# 4 Oligo.3 ENSG000001… TESPA1            2.01             0.0171             0.366              0.0492 TRUE      TRUE
# 5 Oligo.3 ENSG000001… CABLES1           2.19             0.0124             0.382              0.0515 TRUE      TRUE
# 6 Oligo.3 ENSG000000… CALCRL            1.79             0.0238            -0.325              0.0588 FALSE     FALSE
# 7 Oligo.3 ENSG000001… CNDP1            -0.686            0.0466            -0.249              0.0971 TRUE      TRUE

# SpX     gene_id    gene_name vlmf_sn_logFC vlmf_sn_adj.P.Val vlmf_xenium_logFC vlmf_xenium_P.Value dir_match validate
# <fct>   <chr>      <chr>             <dbl>             <dbl>             <dbl>               <dbl> <lgl>     <lgl>   
# 1 L6~SpX9 ENSG00000… ENC1              1.61            0.0238              0.332              0.0320 TRUE      TRUE    
# 2 L6~SpX9 ENSG00000… OPALIN           -0.653           0.0350             -0.261              0.0350 TRUE      TRUE    
# 3 L6~SpX9 ENSG00000… PTPRD            -0.493           0.0245             -0.360              0.0337 TRUE      TRUE    
# 4 L6~SpX9 ENSG00000… IGFBP5            2.02            0.0166             -0.368              0.0359 FALSE     FALSE   
# 5 L6~SpX9 ENSG00000… CPNE4             1.87            0.0350              0.361              0.0289 TRUE      TRUE    
# 6 L6~SpX9 ENSG00000… TESPA1            2.01            0.0171              0.361              0.0360 TRUE      TRUE    
# 7 L6~SpX9 ENSG00000… CNDP1            -0.686           0.0466             -0.276              0.0720 TRUE      TRUE    
# 8 L6~SpX9 ENSG00000… LPAR1            -0.779           0.0313             -0.292              0.0707 TRUE      TRUE    
# 9 L6~SpX9 ENSG00000… SLC17A7           1.47            0.0388              0.201              0.0818 TRUE      TRUE    
# 10 L6~SpX9 ENSG00000… CABLES1           2.19            0.0124              0.441              0.0935 TRUE      TRUE    
# 11 L6~SpX9 ENSG00000… FOS               2.77            0.00170             0.886              0.0765 TRUE      TRUE   

# vlmf_data_tb_xenium |> select(cluster, gene_name, vlmf_xenium_t, vlmf_sn_t) |> mutate(dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0)) |> print(n = 20)

#### Validation SpX barplot ####

validation_spX_barplot <- validation_summary |>
    filter(n_validate > 0) |>
    mutate(cluster = droplevels(cluster)) |>
    ggplot(aes(x = SpX, y = n_validate, fill = SpX)) +
    geom_col() +
    geom_text(aes(label = n_validate)) +
    scale_fill_manual(values = SpX_colors) +
    facet_wrap(~cluster, ncol = 1, strip.position="right") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "None") 

ggsave(validation_spX_barplot, filename = here(plot_dir, sprintf("%s_validation_SpX_barplot.png", opt$datatype)), height = 8)
    
#### Xenium vs. snRNA-seq correlation ####

# sn_xenium_all_pairs <- vlmf_data_tb |>
#     select(SpX, gene_id, gene_name, starts_with("vlmf")) |> 
#     rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_")) |>
#     rename(cluster_xenium = cluster) |>
#     inner_join(sn_DEG_data |> rename(cluster_sn = cluster), relationship = "many-to-many")


sn_xenium_cor <- vlmf_data_tb_xenium |>
    mutate(cluster_xenium = cluster_SpX, cluster_sn = cell_type_anno) |>
    group_by(cluster_sn, cluster_xenium) |>
    filter(n() >= 10) |> 
    summarise(
        tidy(cor.test(vlmf_sn_t, vlmf_xenium_t, method="spearman")),
        n_genes=n(),
        .groups="drop"
    ) |>
    select(cluster_sn, cluster_xenium, n_genes, cor = estimate, p_value = p.value)|> 
    left_join(vlmf_data_tb_xenium |>
                  mutate(cluster_xenium = cluster_SpX, cluster_sn = cell_type_anno) |>
                  filter(vlmf_sn_adj.P.Val < 0.05) |>
                  group_by(cluster_sn, cluster_xenium) |>
                  filter(n() >= 10) |> ## filter for atleast 10 genes after vlmf pval filter
                  summarise(
                      tidy(cor.test(vlmf_sn_t, vlmf_xenium_t, method="spearman")),
                      n_genes_fdr05=n(),
                      .groups="drop"
                  ) |>
                  select(cluster_sn, cluster_xenium, n_genes_fdr05, cor_fdr05 = estimate, p_value_fdr05 = p.value)
    ) |>
    mutate(FDR = p.adjust(p_value, method="BH"),
           FDR_fdr05 = p.adjust(p_value_fdr05, method="BH"),
           )

# sn_xenium_cor |> filter(cluster_sn == "Oligo.3")
sn_xenium_cor |> arrange(-cor_fdr05)

summary(sn_xenium_cor$n_genes_fdr05)
summary(sn_xenium_cor$cor_fdr05)

write_csv(sn_xenium_cor, file = here(data_dir, sprintf("%s_SpX_v_sn_tstat_cor.csv", opt$datatype)))

sn_xenium_cor_v_genes <- sn_xenium_cor |>
    ggplot(aes(x = n_genes, y = cor, color = cluster_xenium)) +
    geom_point() +
    # geom_text_repel(aes(label = cluster_xenium))+
    theme_bw()

ggsave(sn_xenium_cor_v_genes, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_v_genes.png", opt$datatype)), width =8)

if(opt$datatype == "Xenium_cell_type_anno_SpX"){
    
    sn_xenium_cor_v_SpX <- sn_xenium_cor |>
        ggplot(aes(x = cluster_xenium, y = cluster_sn)) +
        geom_point(aes(size = -log10(p_value), color = cor)) +
        scale_colour_gradient2() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
    
    ggsave(sn_xenium_cor_v_SpX, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_v_SpX.png", opt$datatype)), width =8)
    
} else if(opt$datatype == "Xenium_Oligo.3_Astro_SpX"){
    
    sn_xenium_cor_v_SpX <- sn_xenium_cor |>
        separate(cluster_xenium, into = c("cluster_xenium", "SpX"), sep = "(?=[^_]*$)", extra = "merge") |>
        ggplot(aes(x = cluster_xenium, y = SpX)) +
        geom_point(aes(size = -log10(p_value), color = cor)) +
        scale_colour_gradient2() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
    
    ggsave(sn_xenium_cor_v_SpX, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_v_SpX.png", opt$datatype)), width =8)
    
}

## summary with matching clusters

validation_summary_cor <- validation_summary |>
    left_join(sn_xenium_cor |> 
                  select(SpX = cluster_xenium,
                         cell_type_anno = cluster_sn,
                         n_genes, cor, p_value, 
                         n_genes_fdr05, cor_fdr05, p_value_fdr05))

write_csv(validation_summary_cor, file = here(data_dir, "xenium_O3_SpX_DEG_validation_summary.csv"))


# sn_xenium_cor_sig_heatmap <- sn_xenium_cor |>
#     select(SpX = cluster_xenium, cor, cor_fdr05, FDR, FDR_fdr05) |>
#     pivot_longer(!SpX) |>
#     separate(name, into = c("metric", "filter")) |>
#     replace_na(list(filter = "None")) |>
#     pivot_wider(names_from = "metric", values_from = "value") |>
#     mutate(sig=case_when(
#         FDR < 0.001 ~ "***",
#         FDR < 0.01  ~ "**",
#         FDR < 0.05  ~ "*",
#         TRUE ~ ""
#     )) |>
#     ggplot(aes(x = filter, y = SpX)) +
#     geom_tile(aes(fill = cor)) +
#     geom_text(aes(label = sig)) +
#     theme_bw() + 
#     scale_fill_gradient2(
#         low = "#2166AC",
#         mid = "white",
#         high = "#D6604D",
#         midpoint = 0
#     ) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# ggsave(sn_xenium_cor_sig_heatmap, filename = here(plot_dir, "xenium_O3_SpX_v_sn_t_stat_cor_sig_heatmap.png"), width =4)


# ## Xenium vs. Xenium
# xenium_xenium_cor <- vlmf_data_tb |>
#     select(cluster, gene_id, vlmf_t) |>
#     pivot_wider(names_from=cluster, values_from=vlmf_t) |>
#     select(-gene_id) |>
#     cor(use="pairwise.complete.obs", method = "spearman") |>
#     reshape2::melt() |>
#     rename(cluster_xenium1 = Var1, cluster_xenium2 = Var2, cor = value)|>
#     mutate(cluster_xenium1 = factor(cluster_xenium1, levels = cluster_levels),
#            cluster_xenium2 = factor(cluster_xenium2, levels = cluster_levels),
#            cluster_match = cluster_xenium1 == cluster_xenium2)
# 
# xenium_xenium_cor_heatmap <- xenium_xenium_cor |>
#     ggplot(aes(x = cluster_xenium1, y = cluster_xenium2)) +
#     geom_tile(aes(fill = cor)) +
#     geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
#     theme_bw() + 
#     scale_fill_gradient2(
#         low = "#2166AC",
#         mid = "white",
#         high = "#D6604D",
#         midpoint = 0
#     ) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# ggsave(xenium_xenium_cor_heatmap, filename = here(plot_dir, "xenium_v_xenium_t_stat_cor_heatmap.png"), width =8)
# 
# # sn vs. sn
# sn_sn_cor <- sn_DEG_data |>
#     select(cluster, gene_id, vlmf_sn_t) |>
#     pivot_wider(names_from=cluster, values_from=vlmf_sn_t) |>
#     select(-gene_id) |>
#     cor(use="pairwise.complete.obs", method = "spearman") |>
#     reshape2::melt() |>
#     rename(cluster_sn1 = Var1, cluster_sn2 = Var2, cor = value) |>
#     mutate(cluster_sn1 = factor(cluster_sn1, levels = cluster_levels),
#            cluster_sn2 = factor(cluster_sn2, levels = cluster_levels),
#            cluster_match = cluster_sn1 == cluster_sn2)
# 
# sn_sn_cor_heatmap <- sn_sn_cor |>
#     ggplot(aes(x = cluster_sn1, y = cluster_sn2)) +
#     geom_tile(aes(fill = cor)) +
#     geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
#     theme_bw() + 
#     scale_fill_gradient2(
#         low = "#2166AC",
#         mid = "white",
#         high = "#D6604D",
#         midpoint = 0
#     ) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# ggsave(sn_sn_cor_heatmap, filename = here(plot_dir, "sn_v_sn_t_stat_cor_heatmap.png"), width =8)
# 
# 
# diff_cor <- sn_xenium_cor |>
#     inner_join(sn_sn_cor |> 
#                   rename(cluster_sn = cluster_sn1, cluster_xenium = cluster_sn2, sn_sn_cor = cor)) |>
#     mutate(cor_diff = cor - sn_sn_cor)
#     
# cor_diff_heatmap <- diff_cor |>
#     ggplot(aes(x = cluster_sn, y = cluster_xenium)) +
#     geom_tile(aes(fill = cor_diff)) +
#     geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
#     theme_bw() + 
#     scale_fill_gradient2(
#         low = "#2166AC",
#         mid = "white",
#         high = "#D6604D",
#         midpoint = 0
#     ) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# ggsave(cor_diff_heatmap, filename = here(plot_dir, "xenium_v_sn_t_stat_cor_diff_heatmap.png"), width =8)
# 

#### compare t-stats scatter plots ####

## custom to filter by cell type - plot all SpX
compare_stats_scatter <- function(dge_tb, cell_type, stat = "t", mX = "vlmf_sn", mY = "vlmf_xenium", FDR_cut_mX = 0.05, pval_cut_mY = 0.05, model_name){

    ## define vars
    statX <- paste0(mX, "_", stat)
    statY <- paste0(mY, "_", stat)
    fdrX <- paste0(mX, "_adj.P.Val")
    pvalY <- paste0(mY, "_P.Value")
    # fdrY <- paste0(mY, "_adj.P.Val")

    # define colors
    signif_colors <- c("purple", "blue", "red")
    # names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "FDR<", FDR_cut_mY))
    names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "P value<", pval_cut_mY))

    # make scatter plot
    stat_scatter <- dge_tb |>
        filter(cell_type_anno == cell_type) |>
        mutate(DE_class = case_when(!!sym(fdrX) < FDR_cut_mX & !!sym(pvalY) < pval_cut_mY~ "sig_both",
                                    !!sym(fdrX) < FDR_cut_mX ~ paste(mX, "FDR<", FDR_cut_mX),
                                    # !!sym(fdrY) < FDR_cut_mY ~ paste(mY, "FDR<", FDR_cut_mY), ## use pval for xenium data
                                    !!sym(pvalY) < pval_cut_mY ~ paste(mY, "P value<", pval_cut_mY),
                                    TRUE ~ "None")) |>
        ggplot(aes(x = !!sym(statX), y = !!sym(statY), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        labs(title = model_name, subtitle = paste(mX, "vs.", mY)) +
        facet_wrap(~SpX) +
        theme_bw()

    plot_fn = sprintf("%s_%s_stat_scatter_%s_%s-v-%s_%s.png", opt$datatype, stat, model_name, mX, mY, cell_type)
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 10, width = 10)

    # return(t_stat_scatter)
}

# compare_stats_scatter(dge_tb = vlmf_data_tb_xenium, cell_type = "Oligo.3", model_name = "carrier")

vlmf_data_tb_xenium$cell_type_anno <- droplevels(vlmf_data_tb_xenium$cell_type_anno)

map(levels(vlmf_data_tb_xenium$cell_type_anno), ~try(compare_stats_scatter(dge_tb = vlmf_data_tb_xenium, cell_type = .x, stat = "t", model_name = "carrier")))

# #### compare stats cluster vs. cluster
# carrier_data_wide_t <- vlmf_data_tb |>
#     select(gene_id, gene_id, SpX, vlmf_t) |>
#     # count(cluster)
#     pivot_wider(values_from = "vlmf_t", names_from = "SpX")
# 
# ggpair_t_stats <- ggpairs(carrier_data_wide_t, columns = 2:ncol(carrier_data_wide_t))
# 
# ggsave(ggpair_t_stats, filename = here(plot_dir, "xenium_O3_SpX_t_stat_ggpairs.png"), height = 15, width = 15)    

#### tstat heatmap Xenium_cell_type_anno ####

if(opt$datatype == "Xenium_cell_type_anno_SpX"){
    
    ## select sn FDR < 0.05 genes we want to plot
    sn_DEG_data_FDR05 <- vlmf_data_tb_xenium |> 
        filter(vlmf_sn_adj.P.Val < 0.05) |>
        select(cell_type_anno, gene_name ,vlmf_sn_logFC, vlmf_sn_t, vlmf_sn_adj.P.Val) |>
        unique() |>
        arrange(cell_type_anno, vlmf_sn_t)
    
    sn_DEG_data_FDR05 |> count(cell_type_anno)
    
    ## t-stat matrix
    t_stat_SpX_mat <- map(unique(sn_DEG_data_FDR05$cell_type_anno), 
                          ~vlmf_data_tb_xenium |>
                              filter(cell_type_anno == .x, vlmf_sn_adj.P.Val < 0.05) |>
                              select(gene_name, SpX, vlmf_xenium_t) |>
                              pivot_wider(values_from = "vlmf_xenium_t", names_from = "SpX") |>
                              column_to_rownames("gene_name") |>
                              as.matrix())
    
    names(t_stat_SpX_mat) <- unique(sn_DEG_data_FDR05$cell_type_anno)
    
    ## pval/signif matrix
    p_SpX_mat <- map(unique(sn_DEG_data_FDR05$cell_type_anno), 
                     ~vlmf_data_tb_xenium |>
                         filter(cell_type_anno == .x, vlmf_sn_adj.P.Val < 0.05) |>
                         select(gene_name, SpX, vlmf_xenium_P.Value) |>
                         pivot_wider(values_from = "vlmf_xenium_P.Value", names_from = "SpX") |>
                         column_to_rownames("gene_name") |>
                         as.matrix())
    
    names(p_SpX_mat) <- unique(sn_DEG_data_FDR05$cell_type_anno)
    
    SpX_signif_sum <- map(p_SpX_mat, ~colSums(.x < 0.1, na.rm = TRUE))
    
    
    signif_SpX_mat <- map(unique(sn_DEG_data_FDR05$cell_type_anno), 
                          ~vlmf_data_tb_xenium |>
                              filter(cell_type_anno == .x, vlmf_sn_adj.P.Val < 0.05) |>
                              mutate(signif = case_when(signif_xenium & validate ~ "X",
                                                        signif_xenium ~ "*",
                                                        TRUE ~ "")) |>
                              select(gene_name, SpX, signif) |>
                              pivot_wider(values_from = "signif", names_from = "SpX") |>
                              column_to_rownames("gene_name") |>
                              as.matrix())
    
    names(signif_SpX_mat) <- unique(sn_DEG_data_FDR05$cell_type_anno)
    
    
    map(t_stat_SpX_mat, dim)
    
    
    ## cell type vs. SpX & APOE annotation data
    spx_cell_prop <- read.csv(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "cell_v_SpX_prop_long.csv"), row.names = 1)
    spx_APOE <- read.csv(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "Xenium_SpX_APOE_mean_logcount.csv"), row.names = 1) |>
        rownames_to_column("SpX") |> arrange(-APOE_mean)
    
    APOE_order <- spx_APOE$SpX
    
    
    ## plot for each cell type
    map(names(signif_SpX_mat), function(ct){
        
        sn_DEG_data_test_df <- sn_DEG_data_FDR05 |> 
            filter(cell_type_anno == ct) |>
            select(gene_name, vlmf_sn_t) |>
            column_to_rownames("gene_name")
        
        
        # row annotation by sn t-stat
        sn_t_col_fun = circlize::colorRamp2(c(min(c(-0.01, sn_DEG_data_test_df$vlmf_sn_t)),
                                              0, 
                                              max(c(0.01, sn_DEG_data_test_df$vlmf_sn_t))), 
                                            colors = c(APOE_carrier_colors_dark[['E2+']], "white", APOE_carrier_colors_dark[['E4+']]))
        
        sn_t_row_ha <- rowAnnotation(
            df = sn_DEG_data_test_df,
            col = list(vlmf_sn_t = sn_t_col_fun)
        )
        
        SpX_order <- SpX_levels[SpX_levels %in% colnames(signif_SpX_mat[[ct]])]
        
        # col annotation by n validation
        SpX_val <- validation_summary |>
            filter(cluster == ct) |>
            ungroup() |>
            mutate(n_opp = n_signif_both - n_validate) |>
            select(SpX, n_validate, n_opp) |>
            column_to_rownames("SpX")
        
        SpX_val <- SpX_val[SpX_order,]
        
        ha_SpX_val = HeatmapAnnotation(n_validated = anno_barplot(as.matrix(SpX_val)))
        
        # col annotation n_cells
        spx_cell_prop_ct <- spx_cell_prop |>
            filter(cell_type_anno == ct) |>
            select(SpX, n_cell) |>
            left_join(spx_APOE, by = join_by(SpX)) |> ## Add APOE expression
            column_to_rownames("SpX")
        
        spx_cell_prop_ct <- spx_cell_prop_ct[SpX_order,]
        
        col_fun_n_cell = circlize::colorRamp2(c(0, max(spx_cell_prop_ct$n_cell)), c("white", "red"))
        col_fun_APOE = circlize::colorRamp2(c(min(1, spx_cell_prop_ct$APOE_mean), max(spx_cell_prop_ct$APOE_mean)), c("white", "purple"))
        
        ha_SpX_cell <- HeatmapAnnotation(df = spx_cell_prop_ct,
                                         col = list(n_cell = col_fun_n_cell,
                                                    APOE_mean = col_fun_APOE))
        
        # reorder tabs 
        signif_SpX_ct <- signif_SpX_mat[[ct]][rownames(sn_DEG_data_test_df),SpX_order, drop = FALSE]
        t_stat_SpX_ct <- t_stat_SpX_mat[[ct]][rownames(sn_DEG_data_test_df),SpX_order, drop = FALSE]
        
        signif_SpX_ct[is.na(signif_SpX_ct)] <- ""
        
        ## color scale for xenium stats
        xenium_t_col_fun = circlize::colorRamp2(
            c(min(c(-0.01,t_stat_SpX_ct), na.rm = TRUE), 
              0, 
              max(c(0.01, t_stat_SpX_ct), na.rm = TRUE)
            ), 
            colors = c(APOE_carrier_colors_dark[['E2+']], "white", APOE_carrier_colors_dark[['E4+']])
        )
        
        pdf(here(plot_dir, sprintf("xenium_SpX_t_stat_heatmap-%s.pdf", ct)), height = 3 + nrow(t_stat_SpX_ct)/5)
        print(
            Heatmap(t_stat_SpX_ct, 
                      name = "xenium\nt-stat",
                      col = xenium_t_col_fun,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      right_annotation =  sn_t_row_ha,
                      top_annotation = ha_SpX_val,
                      bottom_annotation = ha_SpX_cell,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(signif_SpX_ct[i, j], x, y, gp = gpar(fontsize = 10))
                      },
                      column_title = ct))
        dev.off()
        
    })
    
} else if(opt$datatype == "Xenium_Oligo.3_Astro_SpX"){
    
    #### tstat heatmap Xenium_Oligo.3_Astro_SpX ####
    
    ## select sn FDR < 0.05 genes we want to plot - same for each type
    sn_DEG_data_FDR05 <- vlmf_data_tb_xenium |> 
        filter(vlmf_sn_adj.P.Val < 0.05) |>
        select(cell_type_anno, gene_name ,vlmf_sn_logFC, vlmf_sn_t, vlmf_sn_adj.P.Val) |>
        unique() |>
        arrange(cell_type_anno, vlmf_sn_t)
    
    sn_DEG_data_FDR05 |> count(cell_type_anno)
    
    names(cluster_levels) <- cluster_levels
    
    ## t-stat matrix
    t_stat_SpX_mat <- map(cluster_levels, 
                          ~vlmf_data_tb_xenium |>
                              filter(cluster == .x, vlmf_sn_adj.P.Val < 0.05) |>
                              select(gene_name, SpX, vlmf_xenium_t) |>
                              pivot_wider(values_from = "vlmf_xenium_t", names_from = "SpX") |>
                              column_to_rownames("gene_name") |>
                              as.matrix())
    
    ## pval/signif matrix
    p_SpX_mat <- map(cluster_levels, 
                     ~vlmf_data_tb_xenium |>
                         filter(cluster == .x, vlmf_sn_adj.P.Val < 0.05) |>
                         select(gene_name, SpX, vlmf_xenium_P.Value) |>
                         pivot_wider(values_from = "vlmf_xenium_P.Value", names_from = "SpX") |>
                         column_to_rownames("gene_name") |>
                         as.matrix())
    
    SpX_signif_sum <- map(p_SpX_mat, ~colSums(.x < 0.1, na.rm = TRUE))
    
    
    signif_SpX_mat <- map(cluster_levels, 
                          ~vlmf_data_tb_xenium |>
                              filter(cluster == .x, vlmf_sn_adj.P.Val < 0.05) |>
                              mutate(signif = case_when(signif_xenium & validate ~ "X",
                                                        signif_xenium ~ "*",
                                                        TRUE ~ "")) |>
                              select(gene_name, SpX, signif) |>
                              pivot_wider(values_from = "signif", names_from = "SpX") |>
                              column_to_rownames("gene_name") |>
                              as.matrix())
    
    map(t_stat_SpX_mat, dim)
    
    
    ## cell type vs. SpX & APOE annotation data
    O3_nnA_prop <- read.csv(here("processed-data", "21_Xenium", "20_xenium_Oligo3_Astro", "Oligo3_Astro_neighbor_SpX_summary.csv"), row.names = 1)
    
    spx_APOE <- read.csv(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "Xenium_SpX_APOE_mean_logcount.csv"), row.names = 1) |>
        rownames_to_column("SpX") |> 
        arrange(-APOE_mean)
    
    APOE_order <- spx_APOE$SpX
    
    
    ## plot for each cell type
    map(cluster_levels, function(ct){
        
        sn_DEG_data_test_df <- sn_DEG_data_FDR05 |> 
            filter(gene_name %in% rownames(t_stat_SpX_mat[[ct]])) |>
            select(gene_name, vlmf_sn_t) |>
            column_to_rownames("gene_name")
        
        
        # row annotation by sn t-stat
        sn_t_col_fun = circlize::colorRamp2(c(min(c(-0.01, sn_DEG_data_test_df$vlmf_sn_t)),
                                              0, 
                                              max(c(0.01, sn_DEG_data_test_df$vlmf_sn_t))), 
                                            colors = c(APOE_carrier_colors_dark[['E2+']], "white", APOE_carrier_colors_dark[['E4+']]))
        
        sn_t_row_ha <- rowAnnotation(
            df = sn_DEG_data_test_df,
            col = list(vlmf_sn_t = sn_t_col_fun)
        )
        
        SpX_order <- SpX_levels[SpX_levels %in% colnames(signif_SpX_mat[[ct]])]
        
        # col annotation by n validation
        SpX_val <- validation_summary |>
            filter(cluster == ct) |>
            ungroup() |>
            mutate(n_opp = n_signif_both - n_validate) |>
            select(SpX, n_validate, n_opp) |>
            column_to_rownames("SpX")
        
        SpX_val <- SpX_val[SpX_order,]
        
        ha_SpX_val = HeatmapAnnotation(n_validated = anno_barplot(as.matrix(SpX_val)))
        
        # col annotation n_cells
        spx_cell_prop_ct <- O3_nnA_prop |>
            filter(Oligo.3_Astro == ct) |>
            select(SpX = reference_SpX, n_cell = n) |>
            left_join(spx_APOE, by = join_by(SpX)) |> ## Add APOE expression
            column_to_rownames("SpX")
        
        spx_cell_prop_ct <- spx_cell_prop_ct[SpX_order,]
        
        col_fun_n_cell = circlize::colorRamp2(c(0, max(spx_cell_prop_ct$n_cell)), c("white", "red"))
        col_fun_APOE = circlize::colorRamp2(c(min(spx_cell_prop_ct$APOE_mean), max(spx_cell_prop_ct$APOE_mean)), c("white", "purple"))
        
        ha_SpX_cell <- HeatmapAnnotation(df = spx_cell_prop_ct,
                                         col = list(n_cell = col_fun_n_cell,
                                                    APOE_mean = col_fun_APOE))
        
        # reorder tabs 
        signif_SpX_ct <- signif_SpX_mat[[ct]][rownames(sn_DEG_data_test_df),SpX_order, drop = FALSE]
        t_stat_SpX_ct <- t_stat_SpX_mat[[ct]][rownames(sn_DEG_data_test_df),SpX_order, drop = FALSE]
        
        signif_SpX_ct[is.na(signif_SpX_ct)] <- ""
        
        ## color scale for xenium stats
        xenium_t_col_fun = circlize::colorRamp2(
            c(min(c(-0.01,t_stat_SpX_ct), na.rm = TRUE), 
              0, 
              max(c(0.01, t_stat_SpX_ct), na.rm = TRUE)
            ), 
            colors = c(APOE_carrier_colors_dark[['E2+']], "white", APOE_carrier_colors_dark[['E4+']])
        )
        
        sn_reg <- ifelse(sn_DEG_data_test_df > 0, "upreg", "downreg")
        
        pdf(here(plot_dir, sprintf("%s_t_stat_heatmap-%s.pdf", opt$datatype, ct)), height = 3 + nrow(t_stat_SpX_ct)/5)
        
        print(Heatmap(t_stat_SpX_ct, 
                      name = "xenium\nt-stat",
                      col = xenium_t_col_fun,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      row_split = sn_reg,
                      right_annotation =  sn_t_row_ha,
                      top_annotation = ha_SpX_val,
                      bottom_annotation = ha_SpX_cell,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(signif_SpX_ct[i, j], x, y, gp = gpar(fontsize = 10))
                      },
                      column_title = ct
        )) 
        
        dev.off()
        
    })
    
    
    #### MEGA tstat heatmap Xenium_Oligo.3_Astro_SpX ####
    
    t_stat_SpX_mat_ALL <- vlmf_data_tb_xenium |>
                              filter(vlmf_sn_adj.P.Val < 0.05) |>
                              select(gene_name, cluster_SpX, vlmf_xenium_t) |>
                              pivot_wider(values_from = "vlmf_xenium_t", names_from = "cluster_SpX") |>
                              column_to_rownames("gene_name") |>
                              as.matrix()
    
    sn_DEG_data_test_df <- sn_DEG_data_FDR05 |> 
        filter(gene_name %in% rownames(t_stat_SpX_mat_ALL)) |>
        select(gene_name, vlmf_sn_t) |>
        column_to_rownames("gene_name")

    # row annotation by sn t-stat
    sn_t_col_fun = circlize::colorRamp2(c(min(c(-0.01, sn_DEG_data_test_df$vlmf_sn_t)),
                                          0, 
                                          max(c(0.01, sn_DEG_data_test_df$vlmf_sn_t))), 
                                        colors = c(APOE_carrier_colors_dark[['E2+']], "white", APOE_carrier_colors_dark[['E4+']]))
    
    sn_t_row_ha <- rowAnnotation(
        df = sn_DEG_data_test_df,
        col = list(vlmf_sn_t = sn_t_col_fun)
    )
    
    # col annotation n_cells
    spx_cell_prop_ct <- O3_nnA_prop |>
        mutate(cluster_SpX = paste0(Oligo.3_Astro, "_", gsub("~SpX.*", "", reference_SpX))) |>
        select(cluster_SpX, SpX = reference_SpX, Oligo.3_Astro, n_cell = n) |>
        left_join(spx_APOE, by = join_by(SpX)) |> ## Add APOE expression
        filter(cluster_SpX %in% colnames(t_stat_SpX_mat_ALL)) |>
        column_to_rownames("cluster_SpX") |>
        arrange(SpX)
    
    # spx_cell_prop_ct <- spx_cell_prop_ct[SpX_order,]
    
    col_fun_n_cell = circlize::colorRamp2(c(0, max(spx_cell_prop_ct$n_cell)), c("white", "red"))
    col_fun_APOE = circlize::colorRamp2(c(min(spx_cell_prop_ct$APOE_mean), max(spx_cell_prop_ct$APOE_mean)), c("white", "purple"))
    
    ha_SpX_cell <- HeatmapAnnotation(df = spx_cell_prop_ct,
                                     col = list(n_cell = col_fun_n_cell,
                                                APOE_mean = col_fun_APOE,
                                                SpX = SpX_colors,
                                                Oligo.3_Astro = cluster_colors))
    
    all(rownames(spx_cell_prop_ct) %in% colnames(t_stat_SpX_mat_ALL))

    signif_SpX_ct_ALL <- signif_SpX_mat_ALL[rownames(sn_DEG_data_test_df),rownames(spx_cell_prop_ct), drop = FALSE]
    t_stat_SpX_mat_ALL <- t_stat_SpX_mat_ALL[rownames(sn_DEG_data_test_df),rownames(spx_cell_prop_ct), drop = FALSE]
    
    sn_reg <- ifelse(sn_DEG_data_test_df > 0, "upreg", "downreg")
    
    pdf(here(plot_dir, sprintf("%s_t_stat_heatmap-ALL.pdf", opt$datatype, ct)), width = 10, height = 3 + nrow(t_stat_SpX_ct)/5)
    
    print(Heatmap(t_stat_SpX_mat_ALL, 
                  name = "xenium\nt-stat",
                  col = xenium_t_col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  row_split = sn_reg,
                  right_annotation =  sn_t_row_ha,
                  # top_annotation = ha_SpX_val,
                  bottom_annotation = ha_SpX_cell,
                  # cell_fun = function(j, i, x, y, width, height, fill) {
                  #     grid.text(signif_SpX_ct[i, j], x, y, gp = gpar(fontsize = 10))
                  # },
                  column_title = "ALL Xenium_Oligo.3_Astro_SpX"
    )) 
    
    dev.off()
    
}



# Heatmap(t_stat_SpX_mat, 
#         name = "t-stat", 
#         cluster_columns  = FALSE, 
#         cluster_rows = TRUE, 
#         right_annotation = sn_t_row_ha)
# 
# Heatmap(t_stat_SpX_mat, 
#         name = "t-stat", 
#         cluster_columns  = FALSE, 
#         cluster_rows = FALSE, 
#         right_annotation = sn_t_row_ha)
# dev.off()

#### qvalue ####
## p-value distibutions

xenium_p_val_ct_histo <- vlmf_data_tb |>
    ggplot(aes(x = vlmf_P.Value)) +
    geom_histogram(binwidth = 0.05) +
    facet_wrap(~cluster, ncol = 3) +
    theme_bw()

ggsave(xenium_p_val_ct_histo, filename = here(plot_dir, "xenium_p_val_histogram.png"), width = 10, height = 10)


sn_p_val_ct_histo <- sn_DEG_data |>
    ggplot(aes(x = vlmf_sn_P.Value)) +
    geom_histogram(binwidth = 0.05) +
    facet_wrap(~cluster, ncol = 3) +
    theme_bw()

ggsave(sn_p_val_ct_histo, filename = here(plot_dir, "sn_fine_p_val_histogram.png"), width = 10, height = 10)
    

library(qvalue)

pi1_per_celltype <- vlmf_data_tb_xenium |>
    filter(vlmf_sn_P.Value < 0.1) |>
    group_by(cluster) |>
    filter(n() >= 30) |>
    summarise(
        pi1 = tryCatch(
            1 - pi0est(vlmf_xenium_P.Value, 
                       method = "bootstrap")$pi0,  # bootstrap more robust than smoother
            error = function(e) NA_real_
        ),
        n_discovery = n()
    ) |>
    arrange(desc(pi1))

pi1_per_celltype |>
    left_join(validation_summary_cor) |> 
    select(cluster, pi1, n_discovery, cor, cor_p01, n_signif_both_dir) |>
    mutate(
        replication_tier = case_when(
            pi1 > 0.4  & cor_p01 > 0.3  ~ "Strong",
            pi1 > 0.4  | cor_p01 > 0.3  ~ "Moderate/discordant",
            pi1 > 0.1  | cor_p01 > 0.1  ~ "Weak",
            is.na(pi1)                   ~ "Inconclusive",
            TRUE                         ~ "No replication"
        )
    )

# cluster          pi1 n_discovery
# <fct>          <dbl>       <int>
# 1 Excit.L2_5.2  0.655           42
# 2 Astro.2       0.392           64
# 3 Oligo.4       0.335           46
# 4 Oligo.5       0.248           66
# 5 Vasc.Endo     0.239           32
# 6 Inhib.Vip     0.199           61
# 7 Excit.L2_5.1  0.0760          65
# 8 Oligo.2       0.0136          42
# 9 Micro.2       0               31
# 10 Oligo.3       0              159
# 11 Excit.L5.1    0               30
# 12 Astro.1      NA               62
# 13 Astro.3      NA               32
# 14 Oligo.1      NA               42

# slurmjobs::job_single('01.2_compile_DGE_Xenium', create_shell = TRUE, memory = '5G', command = "Rscript 01.2_compile_DGE_Xenium.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()