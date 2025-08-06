## Louise Huuki-Myers, June 2025
## Create summary plots for dge data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("ggrepel")
library("getopt")
library("ComplexHeatmap")

# Import command-line parameters
# scec <- matrix(
#     c("datatype", "d", "1", "character", "Data type",
#       "contrast", "c", "1", "character", "contrast"),
#     ncol = 5, byrow = TRUE
# )
# opt <- getopt(scec)
# 
# datatype <- opt$datatype
# contrast <- opt$contrast

## test
datatype = "sn_broad"
contrast = "ancestry"

# datatype = "sn_fine"
# datatype = "Visium"

## set plot dir
plot_dir <- here("plots", "13_compile_DGE", "11_summary_plots_contrast", paste0("contrast_", contrast))
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
} else if(datatype == "sn_fine"){
    load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
    cluster_colors <- cell_type_colors$anno
    cluster_levels <- names(cell_type_colors$anno)
}else if(datatype == "Visium"){
    load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
    
}

cluster_levels <- cluster_levels[cluster_levels != "Other"]

## contrast data locations
if(contrast == "ancestry"){
    dge_dir = "05_compile_DGE_ancestry"
} else if(contrast == "Sex"){
    dge_dir = "09_compile_DGE_Sex"
} 


#### load data ####
dge_data <- readRDS(here("processed-data", "13_compile_DGE", dge_dir, datatype,
                         sprintf("DGE_results_%s_%s.Rds", contrast, datatype))) |>
    mutate(cluster = factor(cluster, levels = cluster_levels))

#### combined summary barplot ####
dge_count <- dge_data |>
    group_by(cluster, contrast) |>
    summarize(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              Up = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC > 0),
              Down = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC < 0))

if(datatype == "sn_fine"){
    
    dge_count <- dge_count |> filter(n_FDR05 > 1)
    
}


dge_summary_bar_reg <- dge_count |>
    select(-n_FDR05) |>
    mutate(Down = -1*Down) |>
    pivot_longer(!c(cluster, contrast), names_to = "reg", values_to = "n_genes") |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg)) +
    geom_col() +
    geom_text(aes(label = ifelse(n_genes != 0, abs(n_genes), ""))) +
    facet_wrap(~contrast, ncol = 1, strip.position = "right") +
    scale_fill_manual(values = c(Up = APOE_carrier_colors[["E4+"]],
                                 Down = APOE_carrier_colors[["E2+"]])) +
    theme_bw() +
    labs(y = sprintf("%s\nn DE genes (FDR < 0.05)", datatype)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom")

ggsave(dge_summary_bar_reg, filename = here(plot_dir, sprintf("DGE_%s_%s_summary_bar_reg_combined.png", contrast, datatype)), height = 4, width = 5)

#### select stat scatter plots ####
source(here("code", "13_compile_DGE", "compare_stats_scatter.R"))

if(datatype == "Visium"){
    
    compare_contrast_stats(dge_data |> filter(cluster == "WM.uf~Sp09D07"), 
                           datatype = sprintf("%s_%s_%s",contrast, datatype, "WM.uf"),
                           height = 5, width = 4)
    
} else if(datatype == "sn_broad"){
    
    compare_contrast_stats(dge_data |> filter(cluster == "Astro"), 
                           datatype = sprintf("%s_%s_%s",contrast, datatype, "Astro"),
                           height = 6, width = 5)
    
}else if(datatype == "sn_fine"){
    
    compare_contrast_stats(dge_data |> filter(cluster == "Oligo.3"), 
                           datatype = sprintf("%s_%s_%s",contrast, datatype, "Oligo.3"),
                           height = 6, width = 5)
    
}


#### logFC heatmaps ####
source(here("code", "13_compile_DGE", "logFC_heatmap.R"))

## top DGEs
topDEGs <- dge_data |>
    group_by(cluster,contrast) |>
    filter(vlmf_adj.P.Val < 0.05) |>
    arrange(vlmf_adj.P.Val) |>
    slice(1:10) |>
    pull(gene_name) |>
    unique()

length(topDEGs)

logFC_Heatmap_contrast(dge_data, gene_list = topDEGs, title = sprintf("topDEGs_%s_%s", contrast, datatype))

## Risk gene heatmap
logFC_Heatmap_contrast(dge_data, AD_risk$symbol, title = sprintf("ADrisk_%s_%s", contrast, datatype))

if(datatype == "cell_type_fine"){
    
    ## carrier
    top_oligo_DEGs <- dge_data |>
        filter(vlmf_adj.P.Val < 0.05, cluster == "Oligo.3") |>
        group_by(cluster, vlmf_logFC > 0) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1:25) |>
        pull(gene_name) |>
        unique()
    
    logFC_Heatmap(data = dge_data |>
                      filter(grepl("Oligo", cluster)), 
                  gene_list = top_oligo_DEGs, 
                  title = "topDEGs_Oligo.3", 
                  cluster_col = TRUE)
    
    logFC_Heatmap(data = dge_data |>
                      filter(grepl("Oligo", cluster)), 
                  gene_list = AD_risk$symbol, 
                  title = "ADrisk_Oligo")
    
    ## carrier + anc
    top_oligo_DEGs_anc <- dge_anc_data |>
        filter(vlmf_adj.P.Val < 0.05, cluster == "Oligo.3") |>
        group_by(cluster, contrast) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1:25) |>
        pull(gene_name) |>
        unique()
    
    logFC_Heatmap(data = dge_anc_data |>
                      filter(grepl("Oligo", cluster)) |>
                      mutate(cluster = gsub("carrier_", "", paste(cluster, contrast))), 
                  gene_list = top_oligo_DEGs_anc, 
                  title = "topDEGs_Anc_Oligo.3",
                  cluster_col = TRUE)
    
    logFC_Heatmap(data = dge_data |>
                      filter(grepl("Oligo", cluster)), 
                  gene_list = AD_risk$symbol, 
                  title = "ADrisk_Oligo")
    
    
}


# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "11_summary_plots_contrast",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

