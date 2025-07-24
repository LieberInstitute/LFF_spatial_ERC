## Louise Huuki-Myers, June 2025
## Create figure friendly volcano plots

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

datatype <- opt$datatype
## test
# datatype = "sn_broad"
# datatype = "sn_fine"
# datatype = "Visium"

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "06_Volcano_plots", datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

#### colors and factors ####

cluster_colors <- c()
cluster_levels <- c()

if(datatype == "sn_broad"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
}else if(datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
}

#### load data ####
dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype,
                         sprintf("DGE_results_carrier_%s.Rds", datatype)))
    mutate(cluster = factor(cluster, levels = cluster_levels))

dge_anc_data <- readRDS(here("processed-data", "13_compile_DGE", "05_compile_DGE_ancestry", datatype,
                         sprintf("DGE_results_ancestry_%s.Rds", datatype)))

dge_data_combined <- bind_rows(dge_data, dge_anc_data) |>
    mutate(cluster = factor(cluster, levels = cluster_levels))

#### combined summary barplot ####
dge_count_combined <- dge_data_combined |>
    group_by(cluster, contrast) |>
    summarize(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              Up = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC > 0),
              Down = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC < 0))


dge_summary_bar_reg <- dge_count_combined |>
    select(-n_FDR05) |>
    mutate(Down = -1*Down) |>
    pivot_longer(!c(cluster, contrast), names_to = "reg", values_to = "n_genes") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = ifelse(n_genes != 0, abs(n_genes), ""))) +
    facet_wrap(~contrast, ncol = 1) +
    scale_fill_manual(values = c(Up = APOE_carrier_colors[["E4+"]],
                                 Down = APOE_carrier_colors[["E2+"]])) +
    theme_bw() +
    labs(title = sprintf("DGE - %s", datatype), y = "n DE genes (FDR < 0.05)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "right")

ggsave(dge_summary_bar_reg, filename = here(plot_dir, sprintf("DGE_%s_summary_bar_reg_combined.png", datatype)), height = 5, width = 6)


#### Volcano plots ####
custom_volcano <- function(data, 
                           clus, 
                           FDR_cut = 0.05, 
                           p_col = "vlmf_P.Value", 
                           fdr_col = "vlmf_adj.P.Val", 
                           lfc_col = "vlmf_logFC", 
                           model_name = "carrier",
                           save = TRUE,
                           save_size = 4){
    
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
        theme(legend.position = "right")
    
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

## plot volcanos
walk(unique(dge_data$cluster), ~custom_volcano(dge_data, clus = .x))




