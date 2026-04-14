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
plot_dir <- here("plots", "13_compile_DGE", "08_summary_plots")
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

#### load data ####
dge_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", datatype,
                         sprintf("DGE_results_carrier_%s.Rds", datatype))) |>
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
    geom_text(aes(y = n_genes + (n_genes/abs(n_genes)*40), label = ifelse(n_genes != 0, abs(n_genes), ""), color = reg)) +
    # geom_text(aes(label = ifelse(n_genes != 0, abs(n_genes), ""), color = reg), 
    #                 size = 3) +
    facet_wrap(~contrast, ncol = 1) +
    scale_fill_manual(values = c(Up = APOE_carrier_colors[["E4+"]],
                                 Down = APOE_carrier_colors[["E2+"]])) +    
    scale_color_manual(values = c(Up = APOE_carrier_colors_dark[["E4+"]],
                                 Down = APOE_carrier_colors_dark[["E2+"]])) +
    theme_bw() +
    labs(title = sprintf("DGE - %s", datatype), y = "n DE genes (FDR < 0.05)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "right")

ggsave(dge_summary_bar_reg, filename = here(plot_dir, sprintf("DGE_%s_summary_bar_reg_combined.png", datatype)), height = 5, width = 6)

if(datatype == "sn_fine"){
    ggsave(dge_summary_bar_reg, filename = here(plot_dir, sprintf("DGE_%s_summary_bar_reg_combined_wide.png", datatype)), height = 8, width = 10)
    
}

#### v gene vs. n signif ####


dge_n_v_fdr <- dge_data_combined |> 
    group_by(contrast, cluster) |>
    summarize(n_gene = n(), 
              n_signif = sum(vlmf_adj.P.Val < 0.05))

dge_n_v_fdr_scatter <- dge_n_v_fdr |> 
    ggplot(aes(x = n_gene, y = n_signif, color = cluster)) +
    geom_point() +
    geom_text_repel(aes(label = cluster)) +
    scale_color_manual(values = cluster_colors) +
    facet_wrap(~contrast) +
    scale_y_log10() +
    theme_bw() +
    theme(legend.position = "None")

ggsave(dge_n_v_fdr_scatter, filename = here(plot_dir, sprintf("DGE_%s_gene_v_signif_scatter.png", datatype)), width = 10)

#### logFC heatmaps ####
source(here("code", "13_compile_DGE", "logFC_heatmap.R"))

## top DGEs
topDEGs <- dge_data |>
    group_by(cluster) |>
    filter(vlmf_adj.P.Val < 0.05) |>
    arrange(vlmf_adj.P.Val) |>
    slice(1:5) |>
    pull(gene_name) |>
    unique()

length(topDEGs)

logFC_Heatmap(data = dge_data, gene_list = topDEGs, title = "topDEGs", datatype = datatype)

## Risk gene heatmap
logFC_Heatmap(AD_risk$symbol, title = "ADrisk")

if(datatype == "sn_fine"){
    
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
                  cluster_col = TRUE,
                  datatype = datatype)
    
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
    
    
    #### check Blanchard genes in select cell types ####
    source(here("external-data", "Blanchard2022", "get_Blanchard_gene_list.R"))
    Blanchard_gene_list <- map(Blanchard_gene_list, ~.x[.x %in% dge_data$gene_name])
    map_int(Blanchard_gene_list, length)
    

    map2(Blanchard_gene_list, names(Blanchard_gene_list),  ~logFC_Heatmap(data = dge_data |>
                                                                             filter(grepl("Oligo", cluster) | cluster == "Astro.3"), 
                                                                         gene_list = .x, 
                                                                         title = .y, 
                                                                         # cluster_col = TRUE,
                                                                         order_genes = TRUE,
                                                                         datatype = datatype))
    
    ## annotations + stats 
    
    Blanchard_gene_anno_stats |> filter(gene_name == "PLP1")
    
    Blanchard_gene_anno_stats <- read_csv(here("external-data", "Blanchard2022","Blanchard_gene_anno_stats.csv")) |>
        unique() |>
        arrange(logFC)
    
    Blanchard_gene_anno <- Blanchard_gene_anno_stats |>
        filter(grp2 %in% c("E34_E33_noAD", "E34_E33_AD")) |>
        select(gene_name, grp2, anno, logFC) |>
        pivot_wider(names_from = "grp2", values_from = "logFC") |>
        column_to_rownames("gene_name")
    
    Blanchard_gene_anno |> count(anno)
    
    logfc_col_fun = colorRamp2(c(min(Blanchard_gene_anno_stats$logFC), 0, max(Blanchard_gene_anno_stats$logFC)), 
                               colors = c("blue", "white", "red"))
    
    gene_col_ha<- HeatmapAnnotation(
        df = Blanchard_gene_anno,
        col = list(anno = c(Cholesterol = "red",
                            Myelination = "blue",
                            `Sterol tf` = "yellow"),
                   E34_E33_AD = logfc_col_fun,
                   E34_E33_noAD = logfc_col_fun)
    )
    
    
    logFC_Heatmap(data = dge_data |>
                      filter(grepl("Oligo", cluster)), 
                  gene_list = rownames(Blanchard_gene_anno), 
                  title = "Blancard_anno", 
                  # cluster_col = TRUE,
                  order_genes = FALSE,
                  datatype = datatype,
                  col_anno = gene_col_ha,
                  save = TRUE)
    
    #### check Jakel immune Oligo markers
    
    jakel_Oligo_stats <- readxl::read_xlsx(here("external-data", "Jakel2019","Jakel2019_STab4.xlsx"), sheet = "imOLGs")
    
    jakel_top_imOLG_markers <- jakel_Oligo_stats |> arrange(-avg_logFC) |> slice(1:20) |> pull(gene)
    
    
    dge_data |> 
      filter(cluster == "Oligo.3") |> 
      filter(gene_name %in% jakel_Oligo_stats$gene) |> 
      select(cluster, gene_name, starts_with("vlmf")) |>
      filter(gene_name == "CD74")
    
    
    logFC_Heatmap(data = dge_data |>
                    filter(grepl("Oligo", cluster)), 
                  gene_list = jakel_top_imOLG_markers, 
                  title = "Jakel_imOLG_markers", 
                  # cluster_col = TRUE,
                  order_genes = FALSE,
                  datatype = datatype,
                  # col_anno = gene_col_ha,
                  save = TRUE)
    
    source(here("code", "utils", "custom_volcano.R"))
    
    vol_imOLG <- custom_volcano(dge_data |> filter(gene_name %in% jakel_Oligo_stats$gene), 
                                clus = "Oligo.3", 
                                save = FALSE, 
                                text = TRUE, 
                                highlight_genes = "CD74") +
      geom_point(size = 1) +
      labs(subtitle = "Jakel et al., 2019 imOLG markers")
    
    ggsave(vol_imOLG, 
           filename = here(plot_dir, sprintf("Volcano_%s_%s-%s_imOGL.png", datatype, "carrier", "Oligo.3")), 
           height = 6, width = 6)
    
    immune_genes <- c("STAT1", "IRF1", "B2M", "SERPINA3", "CD74", "HLA-A","HLA-B","HLA-C")
    
    dge_data |> 
      filter(cluster == "Oligo.3") |> 
      filter(gene_name %in% immune_genes) |> 
      select(cluster, gene_name, starts_with("vlmf"))
    
    # cluster gene_name vlmf_logFC vlmf_AveExpr vlmf_t vlmf_P.Value vlmf_adj.P.Val vlmf_B
    # <fct>   <chr>          <dbl>        <dbl>  <dbl>        <dbl>          <dbl>  <dbl>
    # 1 Oligo.3 STAT1          0.803         4.93  3.14       0.00381         0.0464  -2.12
    # 2 Oligo.3 CD74           1.18          3.15  2.35       0.0259          0.121   -3.47
    # 3 Oligo.3 HLA-A          0.818         4.23  1.97       0.0580          0.195   -4.23
    # 4 Oligo.3 HLA-B          0.645         3.76  1.80       0.0815          0.239   -4.48
    # 5 Oligo.3 IRF1           0.732         2.77  1.57       0.127           0.313   -4.65
    # 6 Oligo.3 HLA-C          0.366         3.38  1.18       0.246           0.469   -5.31
    # 7 Oligo.3 B2M            0.177         6.26  0.858      0.398           0.619   -6.32
    
    source(here("code", "13_compile_DGE", "DE_gene_set_enrichment.R"))
    
    imOGL_enrich <- DE_gene_set_enrichment(gene_list = list(imOLG = jakel_Oligo_stats$gene), 
                           de_results = dge_data,
                           fdr_cut = 0.05)
    
    imOGL_enrich |> arrange(-OR) |> head()
    #         OR         Pval NumSig SetSize         test fdr_cut DEG_n    ID  common_genes
    # 1 8.346211 1.940734e-03      4     110   Oligo.5_up    0.05    47 imOLG  ABCC4, PTPRC, SAMSN1, CHST11
    # 2 6.603897 1.247513e-14     33     116   Oligo.3_up    0.05   679 imOLG  ZFP36L1, DOCK8, CAMK1D, MEF2C, SAT1, AKAP13, SRGN, ARHGAP26, PLCB1, CHST11...
}

# 
# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium")),
#                     create_shell = TRUE,
#                     name = "08_summary_plots",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

