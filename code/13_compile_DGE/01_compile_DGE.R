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
# opt$datatype = "Visium"

data_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "01_compile_DGE", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 
main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

#### colors and factors ####

cluster_colors <- c()
cluster_levels <- c()

if(opt$datatype == "sn_broad"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
}else if(opt$datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
}


#### pseudobulkDGE data ####

pseudobulkDGE_fn <- list(sn_broad = here("processed-data", "08_pseudoBulkDGE_sn", "05_pseudoBulkDGE_sn_broad", "sn_pseudoBulkDGE.rds"),
                          Visium = here("processed-data", "09_pseudoBulkDGE_Visium", "03_pseudoBulkDGE_Visium", "Visium_pseudoBulkDGE.rds"))

pseudobulkDGE_fn <- pseudobulkDGE_fn[[opt$datatype]]

stopifnot(file.exists(pseudobulkDGE_fn))

## load main mod data
pseudobulkDGE_data <- readRDS(pseudobulkDGE_fn)[c("carrier", "APOE", "E4E4", "carrier_i", "APOE_i", "E4E4_i")]
## fix names
names(pseudobulkDGE_data) <- main_mods

pseudobulkDGE_data$carrier$Astro

edit_pseudobulkDEG_data <- function(df_list){
    
    cluster_names <- names(df_list)
    
    df_combined <- map2_dfr(df_list, cluster_names, function(df, cluster){
        df2 <- df |>
            as.data.frame() |>
            dplyr::select(gene_name,
                          gene_id,
                          pbDGE_logFC = logFC,
                          pbDGE_AveExpr = AveExpr,
                          pbDGE_t = t,
                          pbDGE_P.Value = P.Value,
                          pbDGE_adj.P.Val = adj.P.Val,
                          pbDGE_B = B) |>
            mutate(cluster = cluster) |>
            as_tibble()
        return(df2)
    })
    return(df_combined)
}

## edit select data
pseudobulkDGE_data_select <- pseudobulkDGE_data[c("carrier", "e4e4")]
pseudobulkDGE_data_tb <- map(pseudobulkDGE_data_select, edit_pseudobulkDEG_data)

map(pseudobulkDGE_data_tb, ~.x |> filter(pbDGE_adj.P.Val < 0.2) |> count(cluster))
map(pseudobulkDGE_data_tb, ~.x |> filter(pbDGE_adj.P.Val < 0.2))

## datatype - mod
pseudobulkDGE_data_tb$carrier
    
#### dreamlet data ####
dreamlet_fn <- list(sn_broad = here("processed-data", "10_dreamlet_sn", "05_compile_dreamlet_sn", "dreamlet_sn_TopTables.Rdata"),
                         Visium = here("processed-data", "11_dreamlet_Visium", "05_compile_dreamlet_Visium", "dreamlet_Visium_TopTables.Rdata"))

dreamlet_fn <- dreamlet_fn[[opt$datatype]]
stopifnot(file.exists(dreamlet_fn))

dreamlet_data <- get(load(dreamlet_fn))

names(dreamlet_data)

head(dreamlet_data$carrier)

if(opt$datatype == "sn_broad"){
    dreamlet_data_tb <- map(dreamlet_data[c("carrier", "e4e4")],
                            ~.x |>
                                as.data.frame() |>
                                as_tibble() |>
                                dplyr::select(gene_name = ID,
                                              cluster = assay,
                                              dream_logFC = logFC,
                                              dream_AveExpr = AveExpr,
                                              dream_t = t,
                                              dream_P.Value = P.Value,
                                              dream_adj.P.Val = adj.P.Val.cell_type,
                                              dream_B = B))
} else if(opt$datatype == "Visium"){
    dreamlet_data_tb <- map(dreamlet_data[c("carrier", "e4e4")],
                                   ~.x |>
                                       as.data.frame() |>
                                       as_tibble() |>
                                       dplyr::select(gene_name = ID,
                                                     cluster = assay,
                                                     dream_logFC = logFC,
                                                     dream_AveExpr = AveExpr,
                                                     dream_t = t,
                                                     dream_P.Value = P.Value,
                                                     dream_adj.P.Val = adj.P.Val.SpD,
                                                     dream_B = B))
}

map(dreamlet_data_tb, colnames)

map(dreamlet_data_tb, ~.x |> filter(dream_adj.P.Val < 0.2) |> count(cluster))
map(dreamlet_data_tb, ~.x |> filter(dream_adj.P.Val < 0.2))

#### voomLmFit data ####

vlmf_dir <-list(sn_broad = "vlmf_sn_broad",
                Visium = "vlmf_Visium")
vlmf_dir <- vlmf_dir[[opt$datatype]]

vlmf_fn <- list.files(here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", vlmf_dir),
                      full.names = TRUE, pattern = ".rds")

names(vlmf_fn) <- map_chr(vlmf_fn, ~gsub("voomLmFit_sn_broad_|voomLmFit_Visium_|.rds", "", basename(.x)))

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
                        mutate(cluster = factor(gsub("_", "~", cluster), levels = cluster_levels))
)


map(vlmf_data_tb, ~.x |> filter(vlmf_adj.P.Val < 0.05) |> count(cluster))
map(vlmf_data_tb, ~.x |> filter(vlmf_adj.P.Val < 0.05))

vlmf_model_summary <- map2_dfr(vlmf_data_tb, names(vlmf_data_tb), 
                               ~.x |> 
                                   mutate(mod = .y) |>
                                   group_by(cluster, mod) |>
                                   summarize(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(vlmf_model_summary_bar, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar.png", opt$datatype)))

vlmf_model_summary_bar_reg <- vlmf_model_summary |>
    select(-n_FDR05) |>
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

## poster version
ggsave(vlmf_model_summary_bar_reg, filename = here(plot_dir, sprintf("%s_vlmf_model_summary_bar_reg.png", opt$datatype)), height = 5, width = 10)


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

# carrier_data <- readRDS(here(data_dir, sprintf("DGE_results_carrier_%s.Rds", opt$datatype)))

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


#### compare stats cluster vs. cluster 

# t-stat
carrier_data_wide_t <- carrier_data |>
    select(gene_id, gene_id, cluster, vlmf_t) |>
    # count(cluster)
    pivot_wider(values_from = "vlmf_t", names_from = "cluster")

ggpair_t_stats <- ggpairs(carrier_data_wide_t, columns = 2:ncol(carrier_data_wide_t), aes(alpha = 0.2))
ggsave(ggpair_t_stats, filename = here(plot_dir, sprintf("%s_t_stat_ggpairs.png", opt$datatype)), height = 10, width = 10)   

## log FC
carrier_data_wide_logFC <- carrier_data |>
    dplyr::select(gene_id, gene_id, cluster, vlmf_logFC) |>
    tidyr::pivot_wider(values_from = "vlmf_logFC", names_from = "cluster")

ggpair_logFC <- ggpairs(carrier_data_wide_logFC, columns = 2:ncol(carrier_data_wide_logFC), aes(alpha = 0.2))
ggsave(ggpair_logFC, filename = here(plot_dir, sprintf("%s_logFC_ggpairs.png", opt$datatype)), height = 10, width = 10)    


# carrier_data |> filter(gene_name == "TBC1D3D") |> select(gene_name, cluster, vlmf_logFC,vlmf_adj.P.Val)

# slurmjobs::job_single('01_compile_DGE_sn_broad', create_shell = TRUE, memory = '5G', command = "Rscript 01_compile_DGE --datatype sn_broad")
# slurmjobs::job_single('01_compile_DGE_Visium', create_shell = TRUE, memory = '5G', command = "Rscript 01_compile_DGE --datatype Visium")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()