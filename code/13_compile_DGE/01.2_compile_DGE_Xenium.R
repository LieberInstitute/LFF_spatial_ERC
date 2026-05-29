## Louise Huuki-Myers, May 2026
## Compile and plot all DGE data from Xenium

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("GGally")
# library("getopt")

opt <- list()
opt$datatype  <- "Xenium"

data_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 
main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

#### colors and factors ####
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
cluster_colors <- cell_type_colors$anno
cluster_levels <- names(cell_type_colors$anno)
cell_type_broad_levels <- names(cell_type_colors$broad)
cell_type_broad_levels <- cell_type_broad_levels[cell_type_broad_levels != "Other"]


#### voomLmFit data ####
vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", "vlmf_Xenium"),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub("voomLmFit_sn_fine_|.rds", "", basename(.x)))

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
    )  |>
    mutate(cell_type_broad = factor(jaffelab::ss(cluster,"\\."), cell_type_broad_levels),
           cluster = factor(cluster, levels = cluster_levels)
    ) |>
    as_tibble()

levels(vlmf_data_tb$cluster)
levels(vlmf_data_tb$cell_type_broad)

vlmf_data_tb |> filter(vlmf_P.Value < 0.05) |> count(cluster) |> print(n = 35)

vlmf_model_summary <- vlmf_data_tb |> 
                                   mutate(mod = "carrier") |>
                                   group_by(cell_type_broad, cluster, mod) |>
                                   summarize(n_genes= n(),
                                             n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
                                             n_pval05 = sum(vlmf_P.Value < 0.05),
                                             nUP = sum(vlmf_P.Value < 0.05 & vlmf_logFC > 0),
                                             nDown = sum(vlmf_P.Value < 0.05 & vlmf_logFC < 0))

vlmf_model_summary |> filter(n_FDR05 > 0)

## n signif bar plots
vlmf_model_summary_bar <- vlmf_model_summary |>
    ggplot(aes(x = cluster, y = n_pval05, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n_pval05), vjust=-.5) +
    scale_fill_manual(values = cluster_colors) +
    facet_grid(mod~cell_type_broad, scales = "free_x", space = "free") +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "p-value < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None")

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar.png", opt$datatype)), height = 12, width = 10)

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_pval05, -n_genes, -n_FDR05) |>
    mutate(nDown = -1*nDown) |>
    pivot_longer(!c(cell_type_broad, cluster, mod), names_to = "reg", values_to = "n_genes") |>
    # filter(startsWith(mod, "apoe") | mod %in% c("E4E4", "carrier")) |>
    filter(mod == "carrier") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = abs(n_genes))) +
    facet_grid(mod~cell_type_broad, scales = "free_x", space = "free") +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar_reg.png", opt$datatype)), width = 10)

## save summary
write.csv(vlmf_model_summary, file = here(data_dir, sprintf("vlmf_model_summary_%s.csv", opt$datatype)))

#### vlmf volcano plots ####

custom_volcano <- function(data, Pval_cut = 0.05, model_name){
    
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
    
    map(cell_type_broad_levels, ~data |> 
            filter(cell_type_broad == .x) |>
            custom_volcano(, model_name = paste0(mod_name,"-", .x)) )
    
}

## plot volcanos
custom_volcano_ct(data = vlmf_data_tb, mod_name = paste0(opt$datatype, "-carrier"))


## filter to risk genes
custom_volcano_ct(data = vlmf_data_tb |> filter(gene_name %in% AD_risk$symbol) , mod_name = paste0(opt$datatype, "-carrier-risk"))

walk(c("carrier"), ~custom_volcano_ct(data = vlmf_data_tb |> filter(gene_name %in% AD_risk$symbol) |> filter(gene_name %in% AD_risk$symbol), 
                                           mod_name = paste0(opt$datatype, "-", .x, "-risk")))

#### Save vlmf data ####

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

vlmf_data_tb_xenium <- vlmf_data_tb |>
    rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_"))|> 
    left_join(sn_DEG_data)

## save data
saveRDS(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_carrier_%s.Rds", opt$datatype)))
write.csv(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_carrier_%s.csv", opt$datatype)), row.names = FALSE)

#### Summarize validation experiments ####


#### Xenium vs. snRNA-seq correlation ####

sn_xenium_cor <- vlmf_data_tb |>
    select(cluster, gene_id, gene_name, starts_with("vlmf")) |> 
    rename_with(~ str_replace(.x, "vlmf_", "vlmf_xenium_"), starts_with("vlmf_")) |>
    rename(cluster_xenium = cluster) |>
    inner_join(sn_DEG_data |> rename(cluster_sn = cluster), relationship = "many-to-many") |>
    group_by(cluster_sn, cluster_xenium) |>
    summarise(
        cor = cor(vlmf_sn_t, vlmf_xenium_t, use = "complete.obs", method = "spearman"),
        n_genes  = n(),
        .groups  = "drop"
    ) |>
    mutate(cluster_match = cluster_sn == cluster_xenium)

sn_xenium_cor |>
    group_by(cluster_sn) |> 
    slice_max(cor)|> 
    arrange(-cluster_match, -cor) |> 
    print(n = 35) 

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

ggsave(sn_xenium_cor_heatmap, filename = here(plot_dir, "xenium_v_sn_t_stat_cor_heatmap.png"), width =8)


sn_xenium_cor_boxplot <- sn_xenium_cor |>
    ggplot(aes(x = cluster_match, y = cor)) +
    geom_boxplot() +
    theme_bw()

ggsave(sn_xenium_cor_boxplot, filename = here(plot_dir, "xenium_v_sn_t_stat_cor_boxplot.png"))


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

ggsave(xenium_xenium_cor_heatmap, filename = here(plot_dir, "xenium_v_xenium_t_stat_cor_heatmap.png"), width =8)

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


#### compare t-stats scatter plots ####

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

map(cell_type_broad_levels, ~comapre_stats_scatter(dge_tb = vlmf_data_tb, ctb = .x, model_name = "carrier"))
map(cell_type_broad_levels, ~comapre_stats_scatter(dge_tb = vlmf_data_tb, ctb = .x, stat = "logFC", model_name = "carrier"))


vlmf_data_tb |> filter(cluster == "Oligo.3") |> filter(vlmf_xenium =)

#### compare stats cluster vs. cluster 
carrier_data_wide_t <- vlmf_data_tb |>
    select(gene_id, gene_id, cluster, vlmf_xenium_t) |>
    # count(cluster)
    pivot_wider(values_from = "vlmf_xenium_t", names_from = "cluster")

ggpair_t_stats <- ggpairs(carrier_data_wide_t, columns = 2:ncol(carrier_data_wide_t))

# ggsave(ggpair_t_stats, filename = "sn_fine_t_stat_ggpairs.png", height = 15, width = 15)    

map(cell_type_broad_levels, function(ctb){
    message(ctb)
    ggpair_t_stats <- ggpairs(carrier_data_wide_t, columns = grep(ctb, colnames(carrier_data_wide_t)), size = 0.5, alpha = 0.5)
    ggsave(ggpair_t_stats, filename = here(plot_dir, sprintf("xenium_t_stat_ggpairs-%s.png", ctb)))
    
})





# slurmjobs::job_single('01.2_compile_DGE_Xenium', create_shell = TRUE, memory = '5G', command = "Rscript 01.2_compile_DGE_Xenium.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()