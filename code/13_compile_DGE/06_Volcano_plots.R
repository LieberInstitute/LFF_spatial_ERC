## Louise Huuki-Myers, June 2025
## Create figure friendly volcano plots

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("getopt")

plot_dir <- here("plots", "13_compile_DGE", "06_Volcano_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### colors and factors ####
# 
# cluster_colors <- c()
# cluster_levels <- c()
# 
# if(opt$datatype == "sn_broad"){
#     load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
#     cluster_colors <- cell_type_colors$broad
#     cluster_levels <- names(cell_type_colors$broad)
# }else if(opt$datatype == "Visium"){
#     load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
#     cluster_colors <- SpD_colors
#     cluster_levels <- names(SpD_colors)
# }


custom_volcano <- function(data, 
                           clus, 
                           FDR_cut = 0.05, 
                           p_col = "vlmf_P.Value", 
                           fdr_col = "vlmf_adj.P.Val", 
                           lfc_col = "vlmf_logFC", 
                           model_name = "carrier",
                           save = TRUE,
                           save_size = 5){
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    names(signif_colors) <- c("both", paste("FDR<", FDR_cut) , "abs(logFC)>1" )
    
    volcano <- data |>
        filter(cluster == clus) |>
        mutate(DE_class = case_when(!!sym(fdr_col) < FDR_cut & abs(vlmf_logFC) > 1  ~ "both",
                                    !!sym(fdr_col) < FDR_cut ~ paste("FDR<", FDR_cut),
                                    abs(!!sym(lfc_col)) > 1 ~ "abs(logFC)>1",
                                    TRUE ~ "None")) |>
        ggplot(aes(x = vlmf_logFC, y = -log10(!!sym(p_col)), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        scale_color_manual(values = signif_colors) +
        labs(x = "log(FC)", y = "-log10(P value)") +
        theme_bw() +
        labs(title = clus) +
        theme(legend.position = "bottom")
    
    if(nrow(data) > 1000){
        volcano <- volcano + 
            # geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) 
            geom_text_repel(aes(label = ifelse(!!sym(fdr_col) < FDR_cut, gene_name, "")), size = 1.5) 
    } else {
        volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
    }
    if(save) ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_%s_%s-%s.png", datatype, model_name, clus)), 
                    height = save_size + 0.5, 
                    width = save_size)
    else return(volcano)    
}

#### Visium ####
datatype = "Visium"

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "06_Volcano_plots", datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data 
dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype,
                         sprintf("DGE_results_carrier_%s.Rds", datatype)))

## print all volcano
walk(unique(dge_data$cluster), ~custom_volcano(dge_data, clus = .x))

#### sn broad ####
datatype = "sn_broad"

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "06_Volcano_plots", datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data 
dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype,
                         sprintf("DGE_results_carrier_%s.Rds", datatype)))

## print all volcano
walk(unique(dge_data$cluster), ~custom_volcano(dge_data, clus = .x))


#### sn fine ####
datatype = "sn_fine"

plot_dir <- here("plots", "13_compile_DGE", "06_Volcano_plots", datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data 
dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype,
                         sprintf("DGE_results_carrier_%s.Rds", datatype)))

## print all volcano
walk(unique(dge_data$cluster), ~custom_volcano(dge_data, clus = .x))
walk(unique(dge_data$cluster), ~custom_volcano(dge_data |> filter(gene_name %in% AD_risk$symbol), clus = .x, model_name = "carrier-risk"))



