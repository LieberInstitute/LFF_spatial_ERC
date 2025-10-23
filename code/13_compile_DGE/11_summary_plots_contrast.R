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
datatype <- opt$datatype
contrast <- opt$contrast

## test
# datatype = "sn_broad"
# contrast = "ancestry"
# contrast = "Sex"

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
    
    contrast_1 <- "carrier_AA"
    contrast_2 <- "carrier_EA"
    
} else if(contrast == "Sex"){
    dge_dir = "09_compile_DGE_Sex"
    
    contrast_1 <- "carrier_F"
    contrast_2 <- "carrier_M"
} 


#### load data ####
dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype, sprintf("DGE_results_carrier_%s.Rds", datatype))) |>
    mutate(ALL_DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "Up",
                                    vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "Down",
                                    TRUE ~ "None"),
           ALL_DE_class_cluster = paste0(cluster, "_", ALL_DE_class))  |>
    select(data_type, cluster, gene_name, ALL_DE_class ,ALL_DE_class_cluster) 

dge_data |> count(ALL_DE_class_cluster)

dge_data_contrast <- readRDS(here("processed-data", "13_compile_DGE", dge_dir, datatype,
                                  sprintf("DGE_results_%s_%s.Rds", contrast, datatype))) |>
    mutate(cluster = factor(cluster, levels = cluster_levels),
           DE_class = case_when(vlmf_logFC > 0 & vlmf_adj.P.Val < 0.05 ~ "Up",
                                vlmf_logFC < 0 & vlmf_adj.P.Val < 0.05 ~ "Down",
                                TRUE ~ "None"),
           DE_class_cluster = paste0(gsub("\\.", "-", cluster), "_",DE_class , "_", gsub("carrier_", "", contrast))
    )

dge_data_contrast |> count(DE_class_cluster)

DE_contrast_summary <- dge_data_contrast |> 
    mutate(DE = ifelse(vlmf_adj.P.Val < 0.05, gsub("carrier_", "", contrast),"None")) |>
    group_by(gene_name, cluster) |> 
    summarise(DE = paste(unique(DE), collapse = ", "),
              n = n()) |>
    ungroup() |>
    count(cluster, DE)

DE_contrast_summary_reg <- dge_data_contrast |> 
    group_by(gene_name, cluster) |> 
    summarise(contrast = paste(sort(contrast), collapse = ", "),
              DE_classes = paste(sort(DE_class_cluster), collapse = ", "),
              n = n()) |>
    ungroup() |>
    count(cluster, contrast, DE_classes)


dge_data_contrast |> 
    filter(DE_class != "None") |> 
    count(contrast, cluster)

dge_data_contrast |> 
    filter(DE_class != "None") |> 
    count(contrast, cluster)

#### combined summary barplot ####
dge_count <- dge_data_contrast |> 
    left_join(dge_data, relationship = "many-to-many") |>
    group_by(cluster, contrast, ALL_DE_class) |>
    mutate(cluster = factor(cluster, levels = cluster_levels)) |>
    summarize(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              Up = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC > 0),
              Down = sum(vlmf_adj.P.Val < 0.05 & vlmf_logFC < 0))


dge_data_contrast |>
    filter(vlmf_adj.P.Val < 0.05) |>
    count(contrast) 

dge_data_contrast |>
    filter(vlmf_adj.P.Val < 0.05) |>
    count(cluster, contrast) 

if(datatype == "sn_fine"){
    
    dge_data_contrast |>
        filter(vlmf_adj.P.Val < 0.05,
               cluster == "Oligo.3") |>
        count(cluster, contrast) 
    
    DE_contrast_summary |> filter(cluster == "Oligo.3")
    # cluster DE           n
    # <fct>   <chr>    <int>
    # 1 Oligo.3 AA, EA      24
    # 2 Oligo.3 AA, None   509
    # 3 Oligo.3 None     10754
    # 4 Oligo.3 None, EA   209
    
    o3_mat <-  matrix(c(509, 24, 10754, 209), byrow = TRUE, ncol = 2)
    fisher.test(o3_mat, alt = "greater")
    
#     Fisher's Exact Test for Count Data
# 
#       data:  o3_mat
#       p-value = 0.9999
#       alternative hypothesis: true odds ratio is greater than 1
#       95 percent confidence interval:
#       0.2838639       Inf
#       sample estimates:
#       odds ratio 
#       0.4122245 
    
    DE_contrast_summary_reg |> filter(cluster == "Oligo.3")
    
    dge_count <- dge_count |> filter(n_FDR05 > 1)
    
}


dge_summary_bar_reg <- dge_count |>
    select(-n_FDR05) |>
    mutate(Down = -1*Down) |>
    pivot_longer(!c(cluster, contrast, ALL_DE_class), names_to = "reg", values_to = "n_genes") |>
    mutate(reg2 = ifelse(ALL_DE_class == reg, paste(reg, "All"), reg)) |>
    ggplot(aes(x = cluster, y = n_genes, fill = reg2)) +
    geom_col() +
    geom_text(aes(label = ifelse(n_genes != 0, abs(n_genes), "")), position = position_stack(vjust = .5), size = 2.5) +
    facet_wrap(~contrast, ncol = 1, strip.position = "right") +
    scale_fill_manual(values = c(Up = APOE_carrier_colors[["E4+"]],
                                `Up All` = APOE_carrier_colors_dark[["E4+"]],
                                 Down = APOE_carrier_colors[["E2+"]],
                                 `Down All` = APOE_carrier_colors_dark[["E2+"]]),
                      name = "Reg") +
    theme_bw() +
    labs(y = sprintf("%s\nn DE genes (FDR < 0.05)", datatype)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom")

ggsave(dge_summary_bar_reg, filename = here(plot_dir, sprintf("DGE_%s_%s_summary_bar_reg_combined.png", contrast, datatype)), height = 4, width = 5)

#### select stat scatter plots ####
source(here("code", "13_compile_DGE", "compare_stats_scatter.R"))

if(datatype == "Visium"){
    
    compare_contrast_stats(dge_data_contrast |> filter(cluster == "WM.uf~Sp09D07"), 
                           datatype = sprintf("%s_%s_%s",contrast, datatype, "WM.uf"),
                           height = 5, width = 4,
                           contrast_1 = contrast_1,
                           contrast_2 = contrast_2)
    
} else if(datatype == "sn_broad"){
    
    if(contrast == "ancestry"){
                compare_contrast_stats(dge_data_contrast |> filter(cluster == "Astro"), 
                               datatype = sprintf("%s_%s_%s",contrast, datatype, "Astro"),
                               height = 6, width = 5,
                               contrast_1 = contrast_1,
                               contrast_2 = contrast_2)
            }else if(contrast == "Sex"){
                compare_contrast_stats(dge_data_contrast |> filter(cluster == "Oligo"), 
                               datatype = sprintf("%s_%s_%s",contrast, datatype, "Oligo"),
                               height = 6, width = 5,
                               contrast_1 = contrast_1,
                               contrast_2 = contrast_2)
            }
    
    
}else if(datatype == "sn_fine"){
    
    compare_contrast_stats(dge_data_contrast |> filter(cluster == "Oligo.3"), 
                           datatype = sprintf("%s_%s_%s",contrast, datatype, "Oligo.3"),
                           height = 6, width = 5,
                           contrast_1 = contrast_1,
                           contrast_2 = contrast_2)
    
    if(contrast == "Sex"){
        
        compare_contrast_stats(dge_data_contrast |> filter(cluster == "Astro.4"), 
                               datatype = sprintf("%s_%s_%s",contrast, datatype, "Astro.4"),
                               height = 6, width = 5,
                               contrast_1 = contrast_1,
                               contrast_2 = contrast_2)
        
        compare_contrast_stats(dge_data_contrast |> filter(cluster == "OPC.5"), 
                               datatype = sprintf("%s_%s_%s",contrast, datatype, "OPC.5"),
                               height = 6, width = 5,
                               contrast_1 = contrast_1,
                               contrast_2 = contrast_2)
            }
    
}

contrast_levels <- unique(dge_data_contrast$contrast)

dge_data_contrast_diff <- dge_data_contrast |>
    select(gene_name, gene_id, cluster, contrast, vlmf_logFC) |> 
    pivot_wider(values_from = vlmf_logFC, names_from = contrast) |>
    mutate(diff = !!sym(contrast_levels[[1]]) - !!sym(contrast_levels[[2]]),
           abs_diff = abs(diff)) 

# dge_data_contrast_diff |> filter(cluster == "Astro") |>arrange(-abs_diff)
# dge_data_contrast_diff |> filter(cluster == "Astro", gene_name == 'JUP')

#### logFC heatmaps ####
source(here("code", "13_compile_DGE", "logFC_heatmap.R"))

cluster_row = FALSE

## top DGEs
if(datatype == "sn_fine"){
    
    topDEGs <- dge_data_contrast |>
        group_by(cluster,contrast) |>
        filter(vlmf_adj.P.Val < 0.05) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1:5) |>
        pull(gene_name) |>
        unique()
    
    length(topDEGs)
    
    logFC_Heatmap_contrast(dge_data_contrast, gene_list = topDEGs, title = sprintf("topDEGs_%s_%s", contrast, datatype), h = 10)
    
    ## Risk gene heatmap
    logFC_Heatmap_contrast(dge_data_contrast, AD_risk$symbol, title = sprintf("ADrisk_%s_%s", contrast, datatype))
    
} else {
    topDEGs <- dge_data_contrast |>
        group_by(cluster,contrast) |>
        filter(vlmf_adj.P.Val < 0.05) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1:10) |>
        pull(gene_name) |>
        unique()
    
    length(topDEGs)
    
    logFC_Heatmap_contrast(dge_data_contrast |> filter(cluster %in% dge_count$cluster), 
                           gene_list = topDEGs, title = sprintf("topDEGs_%s_%s", contrast, datatype),
                           h = 6)
    
    ## Risk gene heatmap
    logFC_Heatmap_contrast(dge_data_contrast, AD_risk$symbol, title = sprintf("ADrisk_%s_%s", contrast, datatype), h = 12)
    
}


if(datatype == "cell_type_fine"){
    
    ## top Oligo 3
    top_oligo_DEGs <- dge_data_contrast |>
        filter(vlmf_adj.P.Val < 0.05, cluster == "Oligo.3") |>
        group_by(cluster,contrast, vlmf_logFC > 0) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1:10) |>
        pull(gene_name) |>
        unique()
    
    logFC_Heatmap_contrast(dge_data_contrast |>
                               filter(grepl("Oligo", cluster)), 
                           gene_list = top_oligo_DEGs, 
                           title = sprintf("topDEGs_%s_%s_Oligo.3", contrast, datatype), 
                           cluster_col = TRUE)    
    
    ## largest diff Oligo 3
    
    contrast_levels <- unique(dge_data_contrast$contrast)
    
    diff_oligo_DEGs <- dge_data_contrast |>
        filter(vlmf_adj.P.Val < 0.05, cluster == "Oligo.3") |>
        select(gene_name,contrast, vlmf_logFC) |> 
        pivot_wider(values_from = vlmf_logFC, names_from = contrast) |>
        mutate(diff = !!sym(contrast_levels[[1]]) - !!sym(contrast_levels[[2]]),
               abs_diff = abs(diff),
               d_reg = case_when(!!sym(contrast_levels[[1]]) > 1  ~ "dd_c1",
               # d_reg = case_when(!!sym(contrast_levels[[1]]) > 1 & !!sym(contrast_levels[[2]]) < -1 ~ "dd_c1",
                                 !!sym(contrast_levels[[2]]) > 1 & !!sym(contrast_levels[[1]]) < -1 ~ "dd_c2",
                                 TRUE ~ NA)) |>
        
        dge_data_contrast |>
        filter(vlmf_adj.P.Val < 0.05, cluster == "Oligo.3") |>
        select(gene_name,contrast, vlmf_logFC) |>
        count(is.na(vlmf_logFC), contrast)
    
    dge_data_contrast |>
        filter(vlmf_adj.P.Val < 0.05, cluster == "Oligo.3") |>
        select(gene_name,contrast, vlmf_logFC) |> 
        pivot_wider(values_from = vlmf_logFC, names_from = contrast) |>
        ggplot(aes(x = carrier_AA, y = carrier_EA)) +
        geom_point()
        
        arrange(-diff)
    
    
    diff_oligo_DEGs |> count(d_reg)   
    diff_oligo_DEGs |> filter(d_reg == "dd_c1") |> arrange(-carrier_EA)
    
    group_by(cluster,contrast, vlmf_logFC > 0) |>
        arrange(vlmf_adj.P.Val) |>
        slice(1:10) |>
        pull(gene_name) |>
        unique()
    
    logFC_Heatmap_contrast(dge_data_contrast |>
                               filter(grepl("Oligo", cluster)), 
                           gene_list = top_oligo_DEGs, 
                           title = sprintf("topDEGs_%s_%s_Oligo.3", contrast, datatype), 
                           cluster_col = TRUE)
    
    ## risk in Oligo sub-types
    logFC_Heatmap_contrast(dge_data_contrast |>
                               filter(grepl("Oligo", cluster)), 
                           gene_list = AD_risk$symbol, 
                           title = sprintf("ADrisk_%s_%s_Oligo", contrast, datatype))
    
}


# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium"),
#                                  contrast = c("Sex", "ancestry")
#                                  ),
#                     create_shell = TRUE,
#                     name = "11_summary_plots_contrast",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

