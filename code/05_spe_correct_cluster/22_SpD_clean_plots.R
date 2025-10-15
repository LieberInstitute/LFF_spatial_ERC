## Louise Huuki-Myers, March 2025
## Create Main fig plots from SPE object

library("spatialLIBD")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("readxl")
library("jaffelab")
library("cowplot")
library("patchwork")
library("DeconvoBuddies")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

plot_dir <- here("plots", "05_spe_correct_cluster", "22_SpD_clean_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
spe
# dim: 30494 122202 

SpD_colors <- metadata(spe)$SpD_colors

pd <- as.data.frame(colData(spe))
colnames(pd)

#### SpD Summary ####
SpD_summary <- pd |>
    group_by(SpD) |>
    summarise(n_spots = n(),
              prop = n_spots/ncol(spe),
              n_donors = length(unique(sample_id)),
              median_sum_umi = median(sum_umi),
              median_sum_gene = median(sum_gene),
              median_expr_chrM_ratio = median(expr_chrM_ratio),
              median_nuclei = median(Nmask_dark_blue))
# |>
#     left_join(enrichment_stats_top_list)


write.csv(SpD_summary, file = here(data_dir, "ERC_Visium_summary_SpD.csv"))

#### Donor Summary ####

sample_summary <- pd |>
    group_by(sample_id, round, Visium_slide)  |>
    summarise(n_spots = n(),
              prop = n_spots/ncol(spe),
              n_donors = length(unique(sample_id)),
              median_sum_umi = median(sum_umi),
              median_sum_gene = median(sum_gene),
              median_expr_chrM_ratio = median(expr_chrM_ratio))

summary(sample_summary)

write.csv(sample_summary, file = here(data_dir, "ERC_Visium_summary_sample.csv"))

#### Quality metrics ####

qc_violin_sum_umi <- ggplot(pd, aes(x = SpD, y = sum_umi, fill = SpD)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = SpD_colors) +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    labs(y = "log10(Sum UMI)") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(qc_violin_sum_umi, filename = here(plot_dir, "ERC_Visium_SpD_QC_violin_sum_umi.png"), width = 7, height =4)

qc_violin_detected <- ggplot(pd, aes(x = SpD, y = sum_gene, fill = SpD)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = SpD_colors) +
    theme_bw() +
    labs(y = "Detected Genes") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(qc_violin_detected, filename = here(plot_dir, "ERC_Visium_SpD_QC_violin_detected.png"), width = 7, height =4)


qc_violin_mito <- ggplot(pd, aes(x = SpD, y = expr_chrM_ratio, fill = SpD)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = SpD_colors) +
    theme_bw() +
    labs(y = "Percent Mito") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(qc_violin_mito, filename = here(plot_dir, "ERC_Visium_QC_SpD_violin_Mito_percent.png"), width = 7, height =4)

#### SpD metrics ####
SpD_info <- read_csv(here("processed-data", "05_spe_correct_cluster", "19_SpD_update_spe", "ERC_spe_cluster_info.csv")) |>
    mutate(SpD = factor(SpD, levels = names(SpD_colors)))

# n nuclei bar plot 

SpD_barplot_n_spots <- SpD_info |>
    ggplot(aes(x = SpD, y =n, fill = SpD)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -.5) +
    scale_fill_manual(values = SpD_colors) +
    theme_bw()  +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "n spots")
    
    
ggsave(SpD_barplot_n_spots, filename = here(plot_dir, "ERC_SpD_barplot_n_spots.png"), height = 5, width = 7)

#### n nuclei ####
summary(pd$CNmask_dark_blue)

pd |>
    group_by(SpD) |>
    summarise(median_nuc = median(CNmask_dark_blue),
              n_0_nuc = sum(CNmask_dark_blue == 0),
              p_0_nuc = n_0_nuc/n(),
              n_100_nuc = sum(CNmask_dark_blue > 100),
              p_100_nuc = n_100_nuc/n(),
              max = max(CNmask_dark_blue))

# SpD           median_nuc n_0_nuc p_0_nuc n_100_nuc p_100_nuc   max
# <fct>              <dbl>   <int>   <dbl>     <int>     <dbl> <int>
# 1 Vasc~Sp09D08           4    1434  0.231        124   0.0199    314
# 2 L1~Sp09D05             3    2134  0.139        243   0.0158    307
# 3 L2.3~Sp09D01           6     801  0.0698        53   0.00462   231
# 4 LD~Sp09D02             3    3871  0.146        244   0.00921   244
# 5 Inhib~Sp09D09          5     316  0.0457        64   0.00926   235
# 6 L5~Sp09D03             6     629  0.0541       104   0.00895   190
# 7 L6~Sp09D04             5    1578  0.0811       281   0.0144    290
# 8 WM.uf~Sp09D07          6    1870  0.149        718   0.0574    398
# 9 WM~Sp09D06             9     419  0.0344       790   0.0649    374

n_nuclei_violin <- ggplot(pd, aes(x = SpD, y = CNmask_dark_blue, fill = SpD)) +
    geom_violin(draw_quantiles = c(.5)) +
    scale_fill_manual(values = SpD_colors) +
    theme_bw() +
    labs(y = "n segmented nuclei") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(n_nuclei_violin, filename = here(plot_dir, "ERC_Visium_SpD_violin_n_nuclei.png"), width = 7, height =4)


#### Spot plots for representative sections ####

table(spe$SpD, spe$sample_id)

rep_sections_tb <- read.csv(here(data_dir, "rep_section.csv")) |>
    filter(rep_section)

single_vis_clus <- vis_clus(
    spe = spe,
    point_size = 1.7,
    colors = SpD_colors,
    sampleid = "Br5517",
    clustervar = "SpD",
    spatial = FALSE,
    guide_point_size = 5
)

ggsave(single_vis_clus, filename = here(plot_dir, "vis_clus_Br5517.png")) 
ggsave(single_vis_clus, filename = here(plot_dir, "vis_clus_Br5517.pdf")) 

cluster_plots <- map(c("AA", "EA"), function(anc) {
# cluster_plots <- map(unique(spe$APOE), function(anc) {
    
    samples <- rep_sections_tb |> filter(Ancestry == anc) |> arrange(APOE)
    # samples <- rep_sections_tb |> filter(APOE == anc) |> arrange(Ancestry)
    
    cluster_row_plots <- map(samples$sample_id, function(s) {
        vis_clus_plot <- vis_clus(
            spe = spe,
            point_size = 1.5,
            colors = SpD_colors,
            sampleid = s,
            clustervar = "SpD"
        ) +
            labs(title = s)  +
            theme(
                legend.position = "None", ## using heat maps label colors
                axis.title.x = element_blank(),
                text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5)
            )
            
        return(vis_clus_plot)
    })
    cluster_row <- Reduce("+", cluster_row_plots) + plot_layout(nrow = 1)
    # ggsave(cluster_row, filename = here(plot_dir, paste0("vis_clust_",k_label,"_row.png")), width = 18)
    return(cluster_row)
})

cluster_grid <- Reduce("/", cluster_plots)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.pdf"), width = 18, height = 9)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.png"), width = 18, height = 9)


#### Plot reduced dims ####
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(spe,
                          prefix = "ERC_spe",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "SpD",
                          color_pal = SpD_colors))


#### Spatial Registration vs. DLPFC & HPC ####
## get reference layer enrichment statistics
layer_modeling_results <- map(c(HumanPilot = "modeling_results", spatialDLPFC = "spatialDLPFC_Visium_modeling_results"), fetch_data)

## Add spatialDLPFC spatial domain annotations to modeling
dlpfc_anno <- read.csv("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_spatial_registration/bayesSpace_layer_annotations.csv") |>
    dplyr::filter(bayesSpace == "k09") |>
    mutate(layer_combo2 = gsub(" ", "~", layer_combo2))

dlpfc_colnames <- colnames(layer_modeling_results$spatialDLPFC$enrichment)
pwalk(dlpfc_anno, function(...) dlpfc_colnames <<- gsub(..5, ..3, dlpfc_colnames))
colnames(layer_modeling_results$spatialDLPFC$enrichment) <- dlpfc_colnames
head(layer_modeling_results$spatialDLPFC$enrichment)

## spatial HPC modeling 
layer_modeling_results$spatialHPC <- list()
layer_modeling_results$spatialHPC$enrichment <- read.csv("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)

## ERC Annotated modeling results
erc_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "modeling_results-SpD.rds"))

colnames(erc_modeling$enrichment) <- gsub("_Sp", "~Sp", colnames(erc_modeling$enrichment))


## correlate 
cor_layer <- map(layer_modeling_results, function(layer_mod){
    
    cor_layer <- layer_stat_cor(stats = erc_modeling$enrichment,
                                modeling_results = layer_mod,
                                model_type = "enrichment",
                                top_n = 100)
    
    cor_layer <- cor_layer[levels(spe$SpD),order(colnames(cor_layer))] ## match factor level for SpD
    
    return(cor_layer)
    
})

anno <- map(cor_layer, ~annotate_registered_clusters(
    cor_stats_layer = .x,
    confidence_threshold = 0.5,
    cutoff_merge_ratio = 0.1
))

anno_summary <- map2_dfr(anno, names(anno), ~.x |> 
                             select(cluster, layer_label) |>
                             mutate(dataset = .y)) |>
    pivot_wider(names_from = "dataset",
                values_from = "layer_label") |>
    arrange(cluster)

write_csv(anno_summary, file = here(data_dir, "ERC_SpD_spatial_registration_anno_summary.csv"))


## save data
save(cor_layer, anno, file = here(data_dir, "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"))

## reference colors
load(here("processed-data","00_project_prep", "spatialDLPFC_Data","spatialDLPFC_SpD_colors.Rdata"), verbose = TRUE)

layer_colors <- list(HumanPilot = spatialLIBD::libd_layer_colors,
                     spatialDLPFC = spatialDLPFC_SpD_colors)
                     
map(names(cor_layer), function(ref){
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s.pdf", ref)), width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = SpD_colors,
        cluster_rows = FALSE,
        cluster_columns = FALSE
    ))
    dev.off()
})

map(names(cor_layer), function(ref){
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_cluster.pdf", ref)), width = 6 + (ncol(cor_layer[[ref]])/7))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = SpD_colors
    ))
    dev.off()
})

#### plot top markers ####

## load pseuodbulked data
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds"))

rownames(spe) <- rowData(spe)$gene_name
rownames(spe_pb) <- rowData(spe_pb)$gene_name

top_DEGs <- sig_genes_extract(n = 10,
                              modeling_results = erc_modeling,
                              model_type = "enrichment",
                              sce_layer = spe_pb) |>
    mutate(cellType.target = test)

plot_marker_express_ALL(
    spe_pb,
    top_DEGs,
    pdf_fn = here(plot_dir, "SpD_top_enrichment_genes_pb.pdf"),
    n_genes = 10,
    rank_col = "top",
    anno_col = NULL,
    gene_col = "gene",
    cellType_col = "SpD",
    color_pal = SpD_colors,
    plot_points = TRUE
)

plot_marker_express_ALL(
    spe,
    top_DEGs,
    pdf_fn = here(plot_dir, "SpD_top_enrichment_genes.pdf"),
    n_genes = 10,
    rank_col = "top",
    anno_col = NULL,
    gene_col = "gene",
    cellType_col = "SpD",
    color_pal = SpD_colors,
    plot_points = FALSE
)


sig_genes_SpD <- readRDS(here("processed-data", "05_spe_correct_cluster", "20_model_pseudobulk_anno", "sig_genes_SpD.rds"))

#### Layer5 sublayer plots ####

sig_genes_SpD |>
    as.data.frame() |>
    dplyr::filter(gene %in% c("ETV1", "BCL11B","FEZF2"),
                  model_type == "enrichment") |>
    arrange(fdr) |>
    group_by(gene) |>
    slice_min(fdr)


L5_sublayer_express <- plot_gene_express(
    spe,
    genes = c("ETV1", "BCL11B","FEZF2"),
    category = "SpD",
    color_pal = SpD_colors,
    ncol = 1
)

ggsave(L5_sublayer_express, filename = here(plot_dir, "L5_sublayer_expression.png"), height = 9, width = 4)

L5_sublayer_express_pb <- plot_gene_express(
    spe_pb,
    genes = c("ETV1", "BCL11B","FEZF2"),
    category = "SpD",
    color_pal = SpD_colors,
    # plot_type = "boxplot",
    plot_points = TRUE,
    ncol = 1
)

ggsave(L5_sublayer_express_pb, filename = here(plot_dir, "L5_sublayer_expression_pb.png"), height = 9, width = 4)

L5_vis_gene <- map(c("ETV1", "BCL11B","FEZF2"), ~vis_gene(spe = spe,
                                                            point_size = 1.7,
                                                            sampleid = "Br5517",
                                                            geneid = .x,
                                                            spatial = TRUE))


ggsave(Reduce("/", L5_vis_gene), filename = here(plot_dir, "L5_sublayer_vis_gene.png"), height = 15, width = 8)


library("escheR")

L5_escheR <- map(c("ETV1", "BCL11B","FEZF2"), function(g){
    
    spe$logcounts <- logcounts(spe)[which(rowData(spe)$gene_name==g),]
    
    p <- make_escheR(spe[, spe$sample_id %in% c("Br5517")]) |>
        add_ground(var = "SpD") |>
        add_fill(var = "logcounts") + 
        scale_fill_gradient(low = "white", high = "black") +
        scale_color_manual(values = SpD_colors) +
        labs(title = g)

}
)

ggsave(Reduce("/", L5_escheR), filename = here(plot_dir, "L5_sublayer_escher.png"), height = 15, width = 8)




# slurmjobs::job_single('22_SpD_clean_plots', create_shell = TRUE, memory = '5G', command = "Rscript 22_SpD_clean_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
