## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("getopt")
library("GGally")

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

data_dir <- here("processed-data", "13_compile_DGE", "05_compile_DGE_ancestry", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "05_compile_DGE_ancestry", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### colors and factors ####

cluster_colors <- c()
cluster_levels <- c()

if(opt$datatype == "sn_broad"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
}else if(opt$datatype == "sn_fine"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$anno
    cluster_levels <- names(cell_type_colors$anno)
    cell_type_broad_levels <- names(cell_type_colors$broad)[names(cell_type_colors$broad) != "Other"]
}else if(opt$datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
}

#### voomLmFit data ####

vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "04_Clusterwise_voomLmFit_ancestry", paste0("vlmf_",opt$datatype)),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub("voomLmFit_ancestry_sn_broad_|voomLmFit_ancestry_sn_fine_|voomLmFit_ancestry_Visium_|.rds", "", basename(.x)))

## read data
vlmf_data <- map(vlmf_fn, readRDS)

## check levels are present
all(names(vlmf_data) %in% cluster_levels)

vlmf_data_tb <- map_dfr(vlmf_data, ~do.call("rbind", .x[c("carrier_AA", "carrier_EA")]) |>
                            dplyr::rename(vlmf_logFC = logFC,
                                          vlmf_AveExpr = AveExpr,
                                          vlmf_t = t,
                                          vlmf_P.Value = P.Value,
                                          vlmf_adj.P.Val = adj.P.Val,
                                          vlmf_B = B
                            )) |>
    mutate(mod = "carrier_by_ancestry") |>
    as_tibble()

if(opt$datatype == "Visium"){
    vlmf_data_tb <- vlmf_data_tb |> mutate(cluster = factor(gsub("_", "~", cluster), levels = cluster_levels))
} else {
    vlmf_data_tb <- vlmf_data_tb |> mutate(cluster = factor(cluster, levels = cluster_levels))
}

levels(vlmf_data_tb$cluster)
stopifnot(!any(is.na(vlmf_data_tb$cluster)))

# vlmf_data_tb|> filter(cluster == "Excit.L2_5.2", gene_name == "MAMLD1") |> select(cluster, gene_name, contrast)
# vlmf_data_tb|> 
#     filter(cluster == "Excit.L5.2", vlmf_adj.P.Val < 0.05) |> 
#     select(cluster, gene_name, contrast, vlmf_adj.P.Val, vlmf_logFC) |>
#     print(n = 175)

vlmf_data_tb|> count(cluster, contrast)
vlmf_data_tb|> filter(vlmf_adj.P.Val < 0.05) |> count(cluster, contrast)
vlmf_data_tb|> arrange(vlmf_adj.P.Val) |> select(cluster, gene_name, vlmf_P.Value, vlmf_adj.P.Val, vlmf_logFC)


(vlmf_model_summary <- vlmf_data_tb |> 
    group_by(cluster, contrast) |>
    summarize(n_genes = n(),
              n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              nUP = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC > 0),
              nDown = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC < 0)))
                               
## n signif bar plots
vlmf_model_summary_bar <- vlmf_model_summary |>
    ggplot(aes(x = cluster, y = n_FDR05, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n_FDR05), vjust=-.5) +
    scale_fill_manual(values = cluster_colors) +
    facet_wrap(~contrast, ncol = 1) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s - carrier by ancestry", opt$datatype), subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None")

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, sprintf("ancestry_%s_vlmf_model_summary_bar.png", opt$datatype)))

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_FDR05) |>
    mutate(nDown = -1*nDown) |>
    pivot_longer(!c(cluster, contrast), names_to = "reg", values_to = "n_genes") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = abs(n_genes))) +
    facet_wrap(~contrast, ncol = 1) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s - carrier by ancestry", opt$datatype), subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, sprintf("ancestry_reg_%s_vlmf_model_summary_bar.png", opt$datatype)))

## save summary
write.csv(vlmf_model_summary, file = here(data_dir, sprintf("vlmf_ancestry_model_summary_%s.csv", opt$datatype)))

#### vlmf volcano plots ####

custom_volcano <- function(data, FDR_cut = 0.05, model_name){
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("both", paste("FDR<", FDR_cut) , "abs(logFC)>1" )
    
    volcano <- data |>
        mutate(DE_class = case_when(vlmf_adj.P.Val < FDR_cut & abs(vlmf_logFC) > 1  ~ "both",
                         vlmf_adj.P.Val < FDR_cut ~ paste("FDR<", FDR_cut),
                         abs(vlmf_logFC) > 1 ~ "abs(logFC)>1",
                         TRUE ~ "None")) |>
        ggplot(aes(x = vlmf_logFC, y = -log10(vlmf_P.Value), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        scale_color_manual(values = signif_colors) +
        facet_grid(contrast~cluster) +
        theme_bw() +
        theme(legend.position = "bottom") +
        labs(title = model_name)
    
    if(nrow(data) > 1000){
        volcano <- volcano + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) 
    } else {
        volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
    }
    ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_plot_%s.png", model_name)), height = 10, width = 12)
}

if(opt$datatype == "sn_fine"){
    
    map(cell_type_broad_levels, ~vlmf_data_tb |>
            filter(grepl(.x, cluster)) |>
            custom_volcano(model_name = paste0("sn_fine-", .x, "-ancestry")))
    
    map(cell_type_broad_levels, ~vlmf_data_tb |> 
            filter(grepl(.x, cluster),
                   gene_name %in% AD_risk$symbol) |>
            custom_volcano(model_name = paste0("sn_fine-", .x, "-ancestry-risk")))
    
} else {
    ## plot volcanos
    custom_volcano(data = vlmf_data_tb, model_name = paste0(opt$datatype, "-ancestry"))
    
    ## filter to risk genes
    custom_volcano(data = vlmf_data_tb|> filter(gene_name %in% AD_risk$symbol), model_name = paste0(opt$datatype, "-ancestry-risk"))
}


#### save data ####

saveRDS(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_ancestry_%s.Rds", opt$datatype)))
write.csv(vlmf_data_tb, file = here(data_dir, sprintf("DGE_results_ancestry_%s.csv", opt$datatype)), row.names = FALSE)

# vlmf_data_tb <- readRDS(here(data_dir, sprintf("DGE_results_ancestry_%s.Rds", opt$datatype)))


#### compare t-stats ####
compare_stats_scatter <- function(dge_tb, stat = "t", mX, mY, FDR_cut_mX = 0.2, FDR_cut_mY = 0.2, model_name){
    
    ## define vars
    statX <- paste0(stat, "_", mX)
    statY <- paste0(stat, "_", mY)
    fdrX <- paste0("fdr_", mX)
    fdrY <- paste0("fdr_", mY)
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "FDR<", FDR_cut_mY))
    
    
    # cor <- dge_tb |>
    #     group_by(cluster) |>
    #     summarise(cor = cor(!!sym(statX), !!sym(statY))) |>
    #     mutate(anno = sprintf("cor=%.2f", cor))
    
    # make scatter plot
    dge_tb_class <- dge_tb |>
        mutate(DE_class = case_when(!!sym(fdrX) < FDR_cut_mX & !!sym(fdrY) < FDR_cut_mY ~ "sig_both",
                                    !!sym(fdrX) < FDR_cut_mX ~ paste(mX, "FDR<", FDR_cut_mX),
                                    !!sym(fdrY) < FDR_cut_mY ~ paste(mY, "FDR<", FDR_cut_mY),
                                    TRUE ~ "None")) 
    
    # return(dge_tb_class)
    
    stat_scatter <- ggplot(data = dge_tb_class, 
                          aes(x = !!sym(statX), y = !!sym(statY), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        # geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        labs(title = model_name, subtitle = paste(mX, "vs.", mY)) + 
        facet_wrap(~cluster) +
        # geom_label(
        #     data = cor, ggplot2::aes(x = -Inf, y = Inf, label = anno),
        #     color = "black",
        #     alpha = 0.5,
        #     vjust = "inward", 
        #     hjust = "inward", 
        #     size = 2.5
        # ) +
        theme_bw() +
        theme(legend.position = "bottom")
    
   return(stat_scatter)
    
}



compare_contrast_stats <- function(dge_tb, stat = "t", m = "vlmf", FDR_cut = 0.05, risk = FALSE, datatype = opt$datatype, text = TRUE, save = TRUE, my_text_size = 1.5){

    ## define vars
    m_stat <- paste0(m, "_", stat)
    fdr <- paste0(m, "_adj.P.Val")
    
    ## create wide data
    
    dge_tb_wide_stat <- dge_tb |>
        select(gene_id, gene_name, cluster, contrast, !!sym(m_stat)) |>
        pivot_wider(names_from = "contrast", values_from = !!sym(m_stat), names_prefix = paste0(stat, "_"))
    
    dge_tb_wide <- dge_tb |>
        select(gene_id, gene_name, cluster, contrast, !!sym(fdr)) |>
        pivot_wider(names_from = "contrast", values_from = !!sym(fdr), names_prefix = "fdr_") |> ## get wide FDR
        full_join(dge_tb_wide_stat, by = join_by(gene_id, gene_name, cluster))
    
    ## filter risk
    if(risk){
        dge_tb_wide <- dge_tb_wide |> filter(gene_name %in% AD_risk$symbol)
    }
    
    ## calc correlation
    cor <- dge_tb_wide |>
        group_by(cluster) |>
        summarise(cor = cor.test(!!sym(paste0(stat,"_carrier_EA")), !!sym(paste0(stat,"_carrier_AA")))$estimate,
                  p_val = cor.test(!!sym(paste0(stat,"_carrier_EA")), !!sym(paste0(stat,"_carrier_AA")))$p.value) |>
        ungroup() |>
        mutate(p_val_adj = p.adjust(p_val, method = "bonf"),
               signif = ifelse(p_val_adj < 0.05, "*", ""),
               anno = sprintf("cor=%.2f%s", cor, signif))
    
    
    ## create scatter plot
    stat_scatter <- compare_stats_scatter(dge_tb = dge_tb_wide, 
                                          stat = stat, 
                                          mX = "carrier_AA", 
                                          mY = "carrier_EA", 
                                          FDR_cut_mX = FDR_cut, 
                                          FDR_cut_mY = FDR_cut, 
                                          model_name = datatype) +
        geom_label(
            data = cor, ggplot2::aes(x = -Inf, y = Inf, label = anno),
            color = "black",
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        )

    if(!text){
        stat_scatter <- stat_scatter 
        plot_fn = sprintf("ancestry_%s_stat_scatter_%s_noText.png", datatype, stat)
    } else if(risk){
        stat_scatter <- stat_scatter + geom_text_repel(aes(label = gene_name), size = my_text_size)
        plot_fn = sprintf("ancestry_%s_stat_scatter_%s_risk.png", datatype, stat)
    } else {
        stat_scatter <- stat_scatter + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = my_text_size)
        plot_fn = sprintf("ancestry_%s_stat_scatter_%s.png", datatype, stat)
    }

    if(save){
        ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 10, width = 10)
        return(cor)
    } else {
        return(stat_scatter)
    }
    
}

compare_carrier_stats <- function(dge_tb, stat = "t", m = "vlmf", FDR_cut = 0.05, risk = FALSE, datatype = opt$datatype, my_text_size = 1.5){
    
    ## load carrier data 
    carrier_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))) |>
        select(cluster, gene_id, t_carrier = vlmf_t, fdr_carrier = vlmf_adj.P.Val)

    ## combine data
    dge_combined <- dge_tb  |>
        select(cluster, contrast, gene_id, gene_name, t_ancestry = vlmf_t, fdr_ancestry = vlmf_adj.P.Val) |>
        left_join(carrier_data, by = join_by(cluster, gene_id))
    
    # dge_combined |> group_by(cluster)|> count(is.na(t_ancestry), is.na(t_carrier))
    # dge_combined |> filter(is.na(t_ancestry) | is.na(t_carrier)) |> group_by(cluster)|> count(is.na(t_ancestry), is.na(t_carrier))
    
    
    ## filter risk
    if(risk){
        dge_combined <- dge_combined |> filter(gene_name %in% AD_risk$symbol)
    }
    
    ## calc correlation
    cor <- dge_combined |> 
        filter(!is.na(t_ancestry) & !is.na(t_carrier)) |>
        group_by(cluster, contrast) |>
        summarise(cor = cor.test(t_ancestry, t_carrier)$estimate,
                  p_val = cor.test(t_ancestry, t_carrier)$p.value) |>
        ungroup() |>
        mutate(p_val_adj = p.adjust(p_val, method = "bonf"),
               signif = ifelse(p_val_adj < 0.05, "*", ""),
               anno = sprintf("cor=%.2f%s", cor, signif))
    
    
    # dge_combined <- dge_combined |> replace_na(list(t_carrier = -Inf, t_ancestry = -Inf))
    
    ## create scatter plot
    stat_scatter <- compare_stats_scatter(dge_tb = dge_combined, 
                                          stat = "t", 
                                          mX = "carrier", 
                                          mY = "ancestry", 
                                          FDR_cut_mX = FDR_cut, 
                                          FDR_cut_mY = FDR_cut, 
                                          model_name = datatype,
                                          my_text_size = my_text_size) +
        facet_grid(cluster~contrast) +
        geom_label(
            data = cor, ggplot2::aes(x = -Inf, y = Inf, label = anno),
            color = "black",
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        )
    
    if(risk){
        stat_scatter <- stat_scatter + geom_text_repel(aes(label = gene_name), size = 1.5)
        plot_fn = sprintf("ancestry_v_carrier_%s_stat_scatter_%s_risk.png", datatype, stat)
    } else {
        stat_scatter <- stat_scatter + geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5)
        plot_fn = sprintf("ancestry_v_carrier_%s_stat_scatter_%s.png", datatype, stat)
    }
    
    
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 12, width = 10)
    
    return(cor)
}

# dge_tb <- vlmf_data_tb |> filter(grepl("Oligo", cluster))

if(opt$datatype == "sn_fine"){

    map(cell_type_broad_levels, ~vlmf_data_tb |>
            filter(grepl(.x, cluster)) |>
            compare_contrast_stats(datatype = paste0(opt$datatype, "_", .x)))
    
    ## compare w/ overall carrier stats
    map(cell_type_broad_levels, ~vlmf_data_tb |>
            filter(grepl(.x, cluster)) |>
            compare_carrier_stats(datatype = paste0(opt$datatype, "_", .x)))    
    
    # 
    # map(cell_type_broad_levels, ~vlmf_data_tb |>
    #         filter(grepl(.x, cluster)) |>
    #         compare_stats_scatter(datatype = paste0(opt$datatype, "_", .x),
    #                               risk = TRUE))
    
    
    ## no text version for Oligo.3
    
    scater_o3_noText <- vlmf_data_tb |>
        filter(cluster == "Oligo.3") |>
        compare_contrast_stats(datatype = paste0(opt$datatype, "_Oligo.3"), text = FALSE, save = FALSE) +
        geom_point(size = 1)
    
    ggsave(scater_o3_noText, filename = here(plot_dir, "ancestry_sn_fine_Oligo.3_stat_scatter_t_noText.png"))
    
    scater_o3_small <- vlmf_data_tb |>
        filter(cluster == "Oligo.3") |>
        compare_contrast_stats(datatype = paste0(opt$datatype, "_Oligo.3"), text = TRUE, save = FALSE, my_text_size = 3) +
        geom_point(size = 1)
    
    ggsave(scater_o3_small, filename = here(plot_dir, "ancestry_sn_fine_Oligo.3_stat_scatter_t_small.png"), height = 7, width = 6)

    
} else {
    # ## compare t-stats
    compare_contrast_stats(vlmf_data_tb)
    compare_contrast_stats(vlmf_data_tb, risk = TRUE)
    
    ## compare w/ overall carrier stats
    compare_carrier_stats(vlmf_data_tb)
    
    ## compare logFC
    # compare_contrast_stats(vlmf_data_tb, stat = "logFC")
}

# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")), 
#                     create_shell = TRUE, 
#                     name = "05_compile_DGE_ancestry", 
#                     create_script = FALSE)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()