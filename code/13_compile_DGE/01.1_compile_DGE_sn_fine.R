## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
# library("getopt")

opt <- list()
opt$datatype  <- "sn_fine"

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


#### voomLmFit data ####
vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", "vlmf_sn_fine"),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub("voomLmFit_sn_fine_|.rds", "", basename(.x)))

## read data
vlmf_data <- map(vlmf_fn, readRDS)

vlmf_data <- list_transpose(vlmf_data)

names(vlmf_data)
head(vlmf_data[[1]][[1]])

vlmf_data_tb <- map(vlmf_data,
                    ~as_tibble(do.call("rbind", .x)) |>
                        dplyr::rename(vlmf_logFC = logFC,
                                      vlmf_AveExpr = AveExpr,
                                      vlmf_t = t,
                                      vlmf_P.Value = P.Value,
                                      vlmf_adj.P.Val = adj.P.Val,
                                      vlmf_B = B
                        )  |>
                        mutate(cell_type_broad = jaffelab::ss(cluster,"\\."),
                               cluster = factor(cluster, levels = cluster_levels)
                               )
)

levels(vlmf_data_tb$carrier$cluster)
unique(vlmf_data_tb$carrier$cell_type_broad)

map(vlmf_data_tb, ~.x |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster))
map(vlmf_data_tb, ~.x |> filter(vlmf_adj.P.Val < 0.05))

vlmf_model_summary <- map2_dfr(vlmf_data_tb, names(vlmf_data_tb), 
                               ~.x |> 
                                   mutate(mod = .y) |>
                                   group_by(cluster, mod) |>
                                   summarize(n_genes= n(),
                                             n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
                                             nUP = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC > 0),
                                             nDown = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC < 0))
                               )

## n signif bar plots
vlmf_model_summary_bar <- vlmf_model_summary |>
    filter(startsWith(mod, "apoe") | mod %in% c("E4E4", "carrier")) |>
    ggplot(aes(x = cluster, y = n_FDR05, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n_FDR05), vjust=-.5) +
    scale_fill_manual(values = cluster_colors) +
    facet_wrap(~mod, ncol = 1) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None")

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar.png", opt$datatype)), height = 12)

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_FDR05, -n_genes) |>
    mutate(nDown = -1*nDown) |>
    pivot_longer(!c(cluster, mod), names_to = "reg", values_to = "n_genes") |>
    # filter(startsWith(mod, "apoe") | mod %in% c("E4E4", "carrier")) |>
    filter(mod == "carrier") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = abs(n_genes))) +
    facet_wrap(~mod, ncol = 1) +
    theme_bw() +
    labs(title = sprintf("voomLmFit - %s", opt$datatype), subtitle = "FDR < 0.05") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar_reg.png", opt$datatype)))

## save summary
write.csv(vlmf_model_summary, file = here(data_dir, sprintf("vlmf_model_summary_%s.csv", opt$datatype)))

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

custom_volcano_ct <- function(data){
    
}

## plot volcanos
walk(c("E4E4", "carrier"), ~custom_volcano(data = vlmf_data_tb[[.x]], model_name = paste0(opt$datatype, "-", .x)))
## filter to risk genes
walk(c("E4E4", "carrier"), ~custom_volcano(data = vlmf_data_tb[[.x]] |> filter(gene_name %in% AD_risk$symbol), 
                                           model_name = paste0(opt$datatype, "-", .x, "-risk")))

#### Add other data to vlmf data ####

carrier_data <- vlmf_data_tb$carrier |>
    left_join(pseudobulkDGE_data_tb$carrier)|>
    left_join(dreamlet_data_tb$carrier)

colnames(carrier_data)

dim(carrier_data)

## save data
saveRDS(carrier_data, file = here(data_dir, sprintf("DGE_results_carrier_%s.Rds", opt$datatype)))
write.csv(carrier_data, file = here(data_dir, sprintf("DGE_results_carrier_%s.csv", opt$datatype)), row.names = FALSE)

#### compare t-stats ####
comapre_stats_scatter <- function(dge_tb, stat = "t", mX, mY, FDR_cut_mX = 0.2, FDR_cut_mY = 0.2, model_name){
    
    ## define vars
    statX <- paste0(mX, "_", stat)
    statY <- paste0(mY, "_", stat)
    fdrX <- paste0(mX, "_adj.P.Val")
    fdrY <- paste0(mY, "_adj.P.Val")

    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "FDR<", FDR_cut_mY))
    
    # make scatter plot
    stat_scatter <- dge_tb |>
        mutate(DE_class = case_when(!!sym(fdrX) < FDR_cut_mX & !!sym(fdrY) < FDR_cut_mY ~ "sig_both",
                                    !!sym(fdrX) < FDR_cut_mX ~ paste(mX, "FDR<", FDR_cut_mX),
                                    !!sym(fdrY) < FDR_cut_mY ~ paste(mY, "FDR<", FDR_cut_mY),
                                    TRUE ~ "None")) |>
        ggplot(aes(x = !!sym(statX), y = !!sym(statY), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        labs(title = model_name, subtitle = paste(mX, "vs.", mY)) + 
        facet_wrap(~cluster) +
        theme_bw()
    
    plot_fn = sprintf("%s_%s_stat_scatter_%s_%s-v-%s.png", opt$datatype, stat, model_name, mX, mY)
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 10, width = 10)
    
    # return(t_stat_scatter)
}

## compare t-stats
comapre_stats_scatter(carrier_data, mX= "dream", mY="pbDGE", model_name = "carrier")
comapre_stats_scatter(carrier_data, mX= "dream", mY="vlmf", model_name = "carrier", FDR_cut_mY = 0.05)
comapre_stats_scatter(carrier_data, mX= "pbDGE", mY="vlmf", model_name = "carrier", FDR_cut_mY = 0.05)

## compare logFC
comapre_stats_scatter(carrier_data, stat = "logFC", mX= "dream", mY="pbDGE", model_name = "carrier")
comapre_stats_scatter(carrier_data, stat = "logFC", mX= "dream", mY="vlmf", model_name = "carrier", FDR_cut_mY = 0.05)
comapre_stats_scatter(carrier_data, stat = "logFC", mX= "pbDGE", mY="vlmf", model_name = "carrier", FDR_cut_mY = 0.05)


#### compare stats cluster vs. clsuter 


# slurmjobs::job_single('01_compile_DGE_sn_broad', create_shell = TRUE, memory = '5G', command = "Rscript 01_compile_DGE --datatype sn_broad")
# slurmjobs::job_single('01_compile_DGE_Visium', create_shell = TRUE, memory = '5G', command = "Rscript 01_compile_DGE --datatype Visium")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()