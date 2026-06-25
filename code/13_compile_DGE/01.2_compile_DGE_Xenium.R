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

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt$datatype  <- "Xenium_Oligo.3_Astro"

data_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 
main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

#### colors and factors ####
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)

if(opt$datatype == "Xenium_cell_type_anno"){

    cluster_levels <- names(cell_type_colors$anno) ## cell type levels
    cluster_colors <- cell_type_colors$anno
    
} else if(opt$datatype == "Xenium_Oligo.3_Astro"){

    cluster_levels <- c("nnA_far_APOE_low", "nnA_far_APOE_high", "nnA_near_APOE_low", "nnA_near_APOE_high")
    
    cluster_colors <- c("nnA_far_APOE_low" = "#87219a",
                        "nnA_far_APOE_high" = "#487800",
                        "nnA_near_APOE_low" = "#b0b3ff",
                        "nnA_near_APOE_high" = "#01d9b4"
                       # "doublet_Oligo3_Astro" = "#ef2396"
                       )
}

cell_type_broad_levels <- names(cell_type_colors$broad)
cell_type_broad_levels <- cell_type_broad_levels[cell_type_broad_levels != "Other"]


#### voomLmFit data ####
message(Sys.time() , "load VLMF data")

vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "06_Clusterwise_voomLmFit_Xenium", paste0("vlmf_", opt$datatype)),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub(sprintf("voomLmFit_%s_|.rds", opt$datatype), "", basename(.x)))

if(opt$datatype == "Xenium_Oligo.3_Astro") vlmf_fn <- vlmf_fn[names(vlmf_fn) %in% cluster_levels] ## exclude doublet data

## read data
vlmf_data <- map(vlmf_fn, readRDS)

# vlmf_data <- list_transpose(vlmf_data)

names(vlmf_data)
head(vlmf_data[[1]])

vlmf_data_tb <- do.call("rbind", vlmf_data) |>
    dplyr::rename(vlmf_logFC = logFC,
                  vlmf_AveExpr = AveExpr,
                  vlmf_t = t,
                  vlmf_P.Value = P.Value,
                  vlmf_adj.P.Val = adj.P.Val,
                  vlmf_B = B
    ) |>
    mutate(cluster = factor(cluster, levels = cluster_levels))  |>
    as_tibble()

levels(vlmf_data_tb$cluster)

if(opt$datatype == "Xenium_cell_type_anno"){
    
    vlmf_data_tb <- vlmf_data_tb |>
        mutate(cell_type_broad = factor(gsub("\\..*", "", cluster), cell_type_broad_levels))
    
    levels(vlmf_data_tb$cell_type_broad)
    
}

vlmf_data_tb |> filter(vlmf_P.Value < 0.10) |> count(cluster) |> print(n = 35)

vlmf_model_summary <- vlmf_data_tb |> 
                                   mutate(mod = "carrier") |>
                                   group_by(cluster, mod) |>
                                   summarize(n_genes= n(),
                                             n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
                                             n_pval10 = sum(vlmf_P.Value < 0.10),
                                             nUP = sum(vlmf_P.Value < 0.05 & vlmf_logFC > 0),
                                             nDown = sum(vlmf_P.Value < 0.05 & vlmf_logFC < 0))

vlmf_model_summary |> filter(n_FDR05 > 0)

## n signif bar plots
vlmf_model_summary_bar <- vlmf_model_summary |>
    ggplot(aes(x = cluster, y = n_pval10, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n_pval10), vjust=-.5) +
    scale_fill_manual(values = cluster_colors) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "p-value < 0.10") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None")

if(opt$datatype == "Xenium_cell_type_anno"){
    vlmf_model_summary_bar <- vlmf_model_summary_bar + facet_grid(mod~cell_type_broad, scales = "free_x", space = "free")
}

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar.png", opt$datatype)), height = 12, width = 10)

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_pval10, -n_genes, -n_FDR05) |>
    mutate(nDown = -1*nDown) |>
    # pivot_longer(!c(cell_type_broad, cluster, mod), names_to = "reg", values_to = "n_genes") |>
    pivot_longer(!c(cluster, mod), names_to = "reg", values_to = "n_genes") |>
    filter(mod == "carrier") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = abs(n_genes))) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

if(opt$datatype == "Xenium_cell_type_anno"){
    vlmf_model_summary_bar_reg <- vlmf_model_summary_bar_reg + facet_grid(mod~cell_type_broad, scales = "free_x", space = "free")
}

ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar_reg.png", opt$datatype)), width = 10)

## save summary
write.csv(vlmf_model_summary, file = here(data_dir, sprintf("vlmf_model_summary_%s.csv", opt$datatype)))

#### vlmf volcano plots ####
message(Sys.time() , "Volcano plots")

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
        facet_wrap(~cluster) +
        theme_bw() +
        labs(title = model_name)
    
    if(nrow(data) > 1000){
        volcano <- volcano + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) 
    } else {
        volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
    }
    ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_plot_%s.png", model_name)), height = 10, width = 10)
}

custom_volcano_ct <- function(data, mod_name){
    
    map(cell_type_broad_levels, ~try(data |> 
            filter(cell_type_broad == .x) |>
            custom_volcano(, model_name = paste0(mod_name,"-", .x)) ))
    
}

if(opt$datatype == "Xenium_cell_type_anno"){
    ## plot volcanos by cell type
    custom_volcano_ct(data = vlmf_data_tb, mod_name = paste0(opt$datatype, "-carrier"))
} else {
    custom_volcano(data = vlmf_data_tb, model_name = paste0(opt$datatype, "-carrier")) 
}


#### Save vlmf data - add snRNA-seq results ####
message(Sys.time() , "Add snRNA-seq stats, Save data")

# carrier_data <- vlmf_data_tb$carrier |>
#     left_join(pseudobulkDGE_data_tb$carrier)|>
#     left_join(dreamlet_data_tb$carrier)
# 
# colnames(carrier_data)
# 
# dim(carrier_data)

## Add sn results
sn_DEG_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds")) |>
    select(cluster, gene_id, gene_name, starts_with("vlmf")) |>
    rename_with(~ str_replace(.x, "vlmf_", "vlmf_sn_"), starts_with("vlmf_"))

if(opt$datatype == "Xenium_cell_type_anno"){
    
    vlmf_data_tb_xenium <- vlmf_data_tb |>
        rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_"))|> 
        left_join(sn_DEG_data) |>
        mutate(signif_xenium = vlmf_xenium_P.Value < 0.1, 
               signif_sn = vlmf_sn_adj.P.Val < 0.05,
               signif_both = signif_xenium & signif_sn,
               dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0),
               validate = signif_both & dir_match)
    
} else if(opt$datatype == "Xenium_Oligo.3_Astro"){
    
    sn_DEG_data <- sn_DEG_data |>
        filter(cluster == "Oligo.3") 
    
    vlmf_data_tb_xenium <- vlmf_data_tb |>
        rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_"))|> 
        left_join(sn_DEG_data|>
                      select(-cluster)) |>
        mutate(signif_xenium = vlmf_xenium_P.Value < 0.1, 
               signif_sn = vlmf_sn_adj.P.Val < 0.05,
               signif_both = signif_xenium & signif_sn,
               dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0),
               validate = signif_both & dir_match)
}


## save data
saveRDS(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_carrier_%s.Rds", opt$datatype)))
write.csv(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_carrier_%s.csv", opt$datatype)), row.names = FALSE)

## save data
saveRDS(vlmf_data_tb_xenium, file = here(data_dir, sprintf("DGE_results_carrier_%s_wSN.Rds", opt$datatype)))
write.csv(vlmf_data_tb_xenium, file = here(data_dir, sprintf("DGE_results_carrier_%s_wSN.csv", opt$datatype)), row.names = FALSE)

#### Summarize validation experiments ####
message(Sys.time() , "Summarize validation experiments")

validation_summary <- vlmf_data_tb_xenium |> 
    filter(!is.na(vlmf_sn_P.Value)) |> 
    group_by(cluster) |>
    summarise(n_signif_xenium = sum(signif_xenium),
              n_signif_sn = sum(signif_sn),
              n_signif_both = sum(signif_both),
              n_validate = sum(signif_both & dir_match),
              percent_valid = 100*sum(validate)/n_signif_sn)

validation_summary |>
    filter(n_signif_sn > 0) |>
    print(n = 32)

if(opt$datatype == "Xenium_cell_type_anno"){
    
    vlmf_data_tb_xenium |> 
        filter(cluster == "Oligo.3", signif_both) |>
        # filter(grepl("Astro", cluster), signif_both) |>
        select(cluster, gene_id, gene_name, vlmf_sn_logFC, vlmf_sn_adj.P.Val, vlmf_xenium_logFC, vlmf_xenium_P.Value, dir_match,validate)
    
    # cluster gene_id     gene_name vlmf_sn_logFC vlmf_sn_adj.P.Val vlmf_xenium_logFC vlmf_xenium_P.Value dir_match validate
    # <fct>   <chr>       <chr>             <dbl>             <dbl>             <dbl>               <dbl> <lgl>     <lgl>   
    # 1 Oligo.3 ENSG000000… CALCRL            1.79             0.0238            -0.443              0.0152 FALSE     FALSE   
    # 2 Oligo.3 ENSG000001… PTPRD            -0.493            0.0245            -0.224              0.0401 TRUE      TRUE    
    # 3 Oligo.3 ENSG000001… CNDP1            -0.686            0.0466            -0.270              0.0466 TRUE      TRUE    
    # 4 Oligo.3 ENSG000001… SLC17A7           1.47             0.0388             0.221              0.0669 TRUE      TRUE    
    # 5 Oligo.3 ENSG000001… TESPA1            2.01             0.0171             0.349              0.0618 TRUE      TRUE    
    # 6 Oligo.3 ENSG000001… STAT4             2.62             0.0330            -0.320              0.0244 FALSE     FALSE   
    # 7 Oligo.3 ENSG000001… SLC24A2          -0.788            0.0224            -0.191              0.0692 TRUE      TRUE    
    # 8 Oligo.3 ENSG000002… NPTXR             1.21             0.0178             0.364              0.0721 TRUE      TRUE    
    # 9 Oligo.3 ENSG000001… ENC1              1.61             0.0238             0.211              0.0900 TRUE      TRUE  
    
    # cluster gene_id     gene_name vlmf_sn_logFC vlmf_sn_adj.P.Val vlmf_xenium_logFC vlmf_xenium_P.Value dir_match validate
    # <fct>   <chr>       <chr>             <dbl>             <dbl>             <dbl>               <dbl> <lgl>     <lgl>   
    # 1 Astro.1 ENSG000002… NPTXR             0.751            0.0310             0.371              0.0149 TRUE      TRUE    
    # 2 Astro.2 ENSG000001… LAMA2            -1.06             0.0364            -0.436              0.0261 TRUE      TRUE    
    # 3 Astro.2 ENSG000001… FZD8              1.01             0.0261             0.337              0.0330 TRUE      TRUE    
    # 4 Astro.2 ENSG000002… TRIL              0.724            0.0440             0.176              0.0772 TRUE      TRUE    
    # 5 Astro.3 ENSG000001… MBP              -1.21             0.0393            -0.400              0.0561 TRUE      TRUE    
    # 6 Astro.3 ENSG000001… SLC24A2          -1.29             0.0126            -0.561              0.0575 TRUE      TRUE
    
    
} else {
    vlmf_data_tb_xenium |> 
        # filter(validate) |>
        filter(gene_name == "OPALIN") |>
        select(cluster, gene_id, gene_name, vlmf_sn_logFC, vlmf_sn_adj.P.Val, vlmf_xenium_logFC, vlmf_xenium_P.Value, vlmf_xenium_t, dir_match,validate) |>
        print(n = 20)
}

# vlmf_data_tb_xenium |> select(cluster, gene_name, vlmf_xenium_t, vlmf_sn_t) |> mutate(dir_match = (vlmf_xenium_t > 0) == (vlmf_sn_t > 0)) |> print(n = 20)

## plot volcanos - just sn signif genes ## some bugs to work out...
# custom_volcano_ct(data = vlmf_data_tb_xenium |> 
#                       filter(signif_sn) |> 
#                       rename(vlmf_P.Value = vlmf_xenium_P.Value,
#                              vlmf_logFC = vlmf_xenium_logFC),
#                   mod_name = paste0(opt$datatype, "-carrier_signif_sn"))


#### Xenium vs. snRNA-seq correlation ####
message(Sys.time() , "Xenium vs. snRNA-seq correlation")

sn_xenium_all_pairs <- vlmf_data_tb |>
    select(cluster, gene_id, gene_name, starts_with("vlmf")) |> 
    rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_")) |>
    rename(cluster_xenium = cluster) |>
    inner_join(sn_DEG_data |> rename(cluster_sn = cluster), relationship = "many-to-many")


sn_xenium_cor <- sn_xenium_all_pairs |>
    group_by(cluster_sn, cluster_xenium) |>
    filter(n() >= 10) |> 
    summarise(
        tidy(cor.test(vlmf_sn_t, vlmf_xenium_t, method="spearman")),
        n_genes=n(),
        .groups="drop"
    ) |>
    select(cluster_sn, cluster_xenium, n_genes, cor = estimate, p_value = p.value)|> 
    left_join(sn_xenium_all_pairs |>
                  filter(vlmf_sn_adj.P.Val < 0.05) |>
                  group_by(cluster_sn, cluster_xenium) |>
                  filter(n() >= 10) |> ## filter for atleast 10 genes after vlmf pval filter
                  summarise(
                      tidy(cor.test(vlmf_sn_t, vlmf_xenium_t, method="spearman")),
                      n_genes_FDR05=n(),
                      .groups="drop"
                  ) |>
                  select(cluster_sn, cluster_xenium, n_genes_FDR05, cor_FDR05 = estimate, p_value_FDR05 = p.value)
    ) |>
    mutate(cluster_match = as.character(cluster_sn) == as.character(cluster_xenium),
           FDR = p.adjust(p_value, method="BH"),
           FDR_fdr05 = p.adjust(p_value_FDR05, method="BH"),
           )

summary(sn_xenium_cor$n_genes_FDR05)

write_csv(sn_xenium_cor, file = here(data_dir, sprintf("%s_v_sn_tstat_cor.csv", opt$datatype)))

sn_xenium_cor |>
    filter(cluster_match) |>
    arrange(cor) |>
    print(n = 31) 
    
sn_xenium_cor |>
    group_by(cluster_sn) |> 
    slice_max(cor)|> 
    arrange(-cluster_match, -cor) |> 
    print(n = 35) 

sn_xenium_cor_v_genes <- sn_xenium_cor |>
    ggplot(aes(x = n_genes, y = cor, color = cluster_match)) +
    geom_point() +
    theme_bw()

ggsave(sn_xenium_cor_v_genes, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_v_genes.png", opt$datatype)), width =8)

## summary with matching clusters

validation_summary_cor <- validation_summary |>
    left_join(sn_xenium_cor |> 
                  filter(cluster_match) |>
                  select(cluster = cluster_sn,
                         n_genes, cor, p_value, 
                         n_genes_FDR05, cor_FDR05, p_value_FDR05))

write_csv(validation_summary_cor, file = here(data_dir, sprintf("%s_DEG_validation_summary.csv", opt$datatype)))


sn_xenium_cor_heatmap <- sn_xenium_cor |>
    ggplot(aes(x = cluster_sn, y = cluster_xenium)) +
    geom_tile(aes(fill = cor)) +
    geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
    theme_bw() + 
    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#D6604D",
        midpoint = 0
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(sn_xenium_cor_heatmap, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_heatmap.png", opt$datatype)), width =8)

sn_xenium_cor_sig_heatmap <- sn_xenium_cor |>
    mutate(sig=case_when(
        FDR < 0.001 ~ "***",
        FDR < 0.01  ~ "**",
        FDR < 0.05  ~ "*",
        TRUE ~ ""
    )) |>
    ggplot(aes(x = cluster_sn, y = cluster_xenium)) +
    geom_tile(aes(fill = cor)) +
    geom_text(aes(label = sig), size = 2) +
    geom_text(aes(label = ifelse(cluster_match, "O", ""))) +
    theme_bw() + 
    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#D6604D",
        midpoint = 0
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(sn_xenium_cor_sig_heatmap, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_sig_heatmap.png", opt$datatype)), width =8)

sn_xenium_cor_FDR05_heatmap <- sn_xenium_cor |>
    ggplot(aes(x = cluster_sn, y = cluster_xenium)) +
    geom_tile(aes(fill = cor_FDR05)) +
    geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
    theme_bw() + 
    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#D6604D",
        midpoint = 0
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(sn_xenium_cor_FDR05_heatmap, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_FDR05_heatmap.png", opt$datatype)), width =8)


sn_xenium_cor_FDR05_sig_heatmap <- sn_xenium_cor |>
    mutate(sig=case_when(
        FDR_fdr05 < 0.001 ~ "***",
        FDR_fdr05 < 0.01  ~ "**",
        FDR_fdr05 < 0.05  ~ "*",
        TRUE ~ ""
    )) |>
    ggplot(aes(x = cluster_sn, y = cluster_xenium)) +
    geom_tile(aes(fill = cor_FDR05)) +
    geom_text(aes(label = sig), size = 2) +
    geom_text(aes(label = ifelse(cluster_match, "O", ""))) +
    theme_bw() + 
    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#D6604D",
        midpoint = 0
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(sn_xenium_cor_FDR05_sig_heatmap, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_FDR05_sig_heatmap.png", opt$datatype)), width =8)


sn_xenium_cor_boxplot <- sn_xenium_cor |>
    ggplot(aes(x = cluster_match, y = cor)) +
    geom_boxplot() +
    theme_bw()

ggsave(sn_xenium_cor_boxplot, filename = here(plot_dir, sprintf("%s_v_sn_t_stat_cor_boxplot.png", opt$datatype)))


## Xenium vs. Xenium
xenium_xenium_cor <- vlmf_data_tb |>
    select(cluster, gene_id, vlmf_t) |>
    pivot_wider(names_from=cluster, values_from=vlmf_t) |>
    select(-gene_id) |>
    cor(use="pairwise.complete.obs", method = "spearman") |>
    reshape2::melt() |>
    rename(cluster_xenium1 = Var1, cluster_xenium2 = Var2, cor = value)|>
    mutate(cluster_xenium1 = factor(cluster_xenium1, levels = cluster_levels),
           cluster_xenium2 = factor(cluster_xenium2, levels = cluster_levels),
           cluster_match = cluster_xenium1 == cluster_xenium2)

xenium_xenium_cor_heatmap <- xenium_xenium_cor |>
    ggplot(aes(x = cluster_xenium1, y = cluster_xenium2)) +
    geom_tile(aes(fill = cor)) +
    geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
    theme_bw() + 
    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#D6604D",
        midpoint = 0
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(xenium_xenium_cor_heatmap, filename = here(plot_dir, sprintf("%s_v_xenium_t_stat_cor_heatmap.png", opt$datatype)), width =8)

if(opt$datatype == "Xenium_cell_type_anno"){
    
    # sn vs. sn
    sn_sn_cor <- sn_DEG_data |>
        select(cluster, gene_id, vlmf_sn_t) |>
        pivot_wider(names_from=cluster, values_from=vlmf_sn_t) |>
        select(-gene_id) |>
        cor(use="pairwise.complete.obs", method = "spearman") |>
        reshape2::melt() |>
        rename(cluster_sn1 = Var1, cluster_sn2 = Var2, cor = value) |>
        mutate(cluster_sn1 = factor(cluster_sn1, levels = cluster_levels),
               cluster_sn2 = factor(cluster_sn2, levels = cluster_levels),
               cluster_match = cluster_sn1 == cluster_sn2)
    
    sn_sn_cor_heatmap <- sn_sn_cor |>
        ggplot(aes(x = cluster_sn1, y = cluster_sn2)) +
        geom_tile(aes(fill = cor)) +
        geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
        theme_bw() + 
        scale_fill_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#D6604D",
            midpoint = 0
        ) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(sn_sn_cor_heatmap, filename = here(plot_dir, "sn_v_sn_t_stat_cor_heatmap.png"), width =8)
    
    
    diff_cor <- sn_xenium_cor |>
        inner_join(sn_sn_cor |> 
                       rename(cluster_sn = cluster_sn1, cluster_xenium = cluster_sn2, sn_sn_cor = cor)) |>
        mutate(cor_diff = cor - sn_sn_cor)
    
    cor_diff_heatmap <- diff_cor |>
        ggplot(aes(x = cluster_sn, y = cluster_xenium)) +
        geom_tile(aes(fill = cor_diff)) +
        geom_text(aes(label = ifelse(cluster_match, "*", ""))) +
        theme_bw() + 
        scale_fill_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#D6604D",
            midpoint = 0
        ) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(cor_diff_heatmap, filename = here(plot_dir, "xenium_v_sn_t_stat_cor_diff_heatmap.png"), width =8)
    
    
}


#### compare t-stats scatter plots ####
message(Sys.time() , "t-stat scatter plots")

comapre_stats_scatter <- function(dge_tb, ctb, stat = "t", mX = "vlmf_sn", mY = "vlmf_xenium", FDR_cut_mX = 0.05, pval_cut_mY = 0.1, model_name){

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
        filter(cell_type_broad == ctb) |>
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
        facet_wrap(~cluster) +
        theme_bw()

    plot_fn = sprintf("%s_%s_stat_scatter_%s_%s-v-%s_%s.png", opt$datatype, stat, model_name, mX, mY, ctb)
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 10, width = 10)

    # return(t_stat_scatter)
}

if(opt$datatype == "Xenium_cell_type_anno"){
    
    map(cell_type_broad_levels, ~comapre_stats_scatter(dge_tb = vlmf_data_tb_xenium, ctb = .x, model_name = "carrier"))
    map(cell_type_broad_levels, ~comapre_stats_scatter(dge_tb = vlmf_data_tb_xenium, ctb = .x, stat = "logFC", model_name = "carrier"))
} else {
    vlmf_data_tb_xenium$cell_type_broad <- "Xenium_Oligo.3_Astro"
    comapre_stats_scatter(dge_tb = vlmf_data_tb_xenium, ctb = "Xenium_Oligo.3_Astro", model_name = "carrier")
}


# vlmf_data_tb |> filter(cluster == "Oligo.3") |> filter(vlmf_xenium =)

#### compare stats cluster vs. cluster 
carrier_data_wide_t <- vlmf_data_tb |>
    select(gene_id, gene_id, cluster, vlmf_t) |>
    # count(cluster)
    pivot_wider(values_from = "vlmf_t", names_from = "cluster")

ggpair_t_stats <- ggpairs(carrier_data_wide_t, columns = 2:ncol(carrier_data_wide_t))

ggsave(ggpair_t_stats, filename = here(plot_dir, sprintf("%s_fine_t_stat_ggpairs.png", opt$datatype)), height = 10, width = 10)

map(cell_type_broad_levels, function(ctb){
    message(ctb)
    ggpair_t_stats <- ggpairs(carrier_data_wide_t, columns = grep(ctb, colnames(carrier_data_wide_t)), size = 0.5, alpha = 0.5)
    ggsave(ggpair_t_stats, filename = here(plot_dir, sprintf("%s_t_stat_ggpairs-%s.png",opt$datatype, ctb)))
    
})


#### qvalue ####
message(Sys.time() , "Check q values")

## p-value distibutions

xenium_p_val_ct_histo <- vlmf_data_tb |>
    ggplot(aes(x = vlmf_P.Value)) +
    geom_histogram(binwidth = 0.05) +
    facet_wrap(~cluster, ncol = 3) +
    theme_bw()

ggsave(xenium_p_val_ct_histo, filename = here(plot_dir, sprintf("%s_p_val_histogram.png", opt$datatype)), width = 10, height = 10)


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
    select(cluster, pi1, n_discovery, cor, cor_FDR05, n_signif_both_dir) |>
    mutate(
        replication_tier = case_when(
            pi1 > 0.4  & cor_FDR05 > 0.3  ~ "Strong",
            pi1 > 0.4  | cor_FDR05 > 0.3  ~ "Moderate/discordant",
            pi1 > 0.1  | cor_FDR05 > 0.1  ~ "Weak",
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


#### t-stat heatmaps ####
message(Sys.time() , "t-stat heatmaps")

if(opt$datatype == "Xenium_Oligo.3_Astro"){
    
    library("ComplexHeatmap")
    
    sn_DEG_data_FDR05 <- vlmf_data_tb_xenium |> 
        filter(vlmf_sn_adj.P.Val < 0.05) |>
        select(gene_name ,vlmf_sn_logFC, vlmf_sn_t, vlmf_sn_adj.P.Val) |>
        unique() |>
        arrange(vlmf_sn_t)
    
    t_stat_SpX_mat <- vlmf_data_tb_xenium |>
                              filter(vlmf_sn_adj.P.Val < 0.05) |>
                              select(gene_name, cluster, vlmf_xenium_t) |>
                              pivot_wider(values_from = "vlmf_xenium_t", names_from = "cluster") |>
                              column_to_rownames("gene_name") |>
                              as.matrix()
    
    p_SpX_mat <- vlmf_data_tb_xenium |>
                         filter(vlmf_sn_adj.P.Val < 0.05) |>
                         select(gene_name, cluster, vlmf_xenium_P.Value) |>
                         pivot_wider(values_from = "vlmf_xenium_P.Value", names_from = "cluster") |>
                         column_to_rownames("gene_name") |>
                         as.matrix()
    
    signif_SpX_mat <- vlmf_data_tb_xenium |>
                              filter(vlmf_sn_adj.P.Val < 0.05) |>
                              mutate(signif = case_when(signif_xenium & validate ~ "X",
                                                        signif_xenium ~ "*",
                                                        TRUE ~ "")) |>
                              select(gene_name, cluster, signif) |>
                              pivot_wider(values_from = "signif", names_from = "cluster") |>
                              column_to_rownames("gene_name") |>
                              as.matrix()
    
    sn_DEG_data_test_df <- sn_DEG_data_FDR05 |> 
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
    
    # col annotation by n validation
    cluster_val <- validation_summary |>
        ungroup() |>
        mutate(n_opp = n_signif_both - n_validate) |>
        select(cluster, n_validate, n_opp) |>
        column_to_rownames("cluster")
    
    cluster_val <- cluster_val[cluster_levels,]
    
    ha_cluster_val = HeatmapAnnotation(n_validated = anno_barplot(as.matrix(cluster_val)))
    
    # reorder tabs 
    signif_SpX_ct <- signif_SpX_mat[rownames(sn_DEG_data_test_df),cluster_levels, drop = FALSE]
    t_stat_SpX_ct <- t_stat_SpX_mat[rownames(sn_DEG_data_test_df),cluster_levels, drop = FALSE]
    
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
    
    pdf(here(plot_dir, sprintf("%s_t_stat_heatmap.pdf", opt$datatype)), height = 3 + nrow(t_stat_SpX_ct)/5)
    print(Heatmap(t_stat_SpX_ct, 
                  name = "xenium\nt-stat",
                  col = xenium_t_col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = TRUE,
                  right_annotation =  sn_t_row_ha,
                  row_split = sn_reg,
                  top_annotation = ha_cluster_val,
                  # bottom_annotation = ha_SpX_cell,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(signif_SpX_ct[i, j], x, y, gp = gpar(fontsize = 10))
                  },
                  column_title = opt$datatype
    ))    
    
    print(Heatmap(t_stat_SpX_ct, 
                  name = "xenium\nt-stat",
                  col = xenium_t_col_fun,
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  right_annotation =  sn_t_row_ha,
                  row_split = sn_reg,
                  top_annotation = ha_cluster_val,
                  # bottom_annotation = ha_SpX_cell,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(signif_SpX_ct[i, j], x, y, gp = gpar(fontsize = 10))
                  },
                  column_title = opt$datatype
    ))
    dev.off()
   
    
}


# slurmjobs::job_single('01.2_compile_DGE_Xenium_cell_type_anno', create_shell = TRUE, memory = '5G', command = "Rscript 01.2_compile_DGE_Xenium.R --datatype Xenium_cell_type_anno")
# slurmjobs::job_single('01.2_compile_DGE_Xenium_Oligo.3_Astro', create_shell = TRUE, memory = '5G', command = "Rscript 01.2_compile_DGE_Xenium.R --datatype Xenium_Oligo.3_Astro")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()