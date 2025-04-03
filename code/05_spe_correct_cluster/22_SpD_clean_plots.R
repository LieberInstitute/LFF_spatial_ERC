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

table(spe$SpD, spe$sample_id)

rep_sections_tb <- read.csv(here(data_dir, "rep_section.csv")) |>
    filter(rep_section)

#### load SpD annotations ####

#### Define colors for SpD ####
load(here("processed-data", "05_spe_correct_cluster", "SpD_colors.Rdata"), verbose = TRUE)

# SpD_colors <- c("Vasc~Sp09D08" = "#E05AD2",
#                 "L1~Sp09D05" = "#0220DE",
#                 "L2.3~Sp09D01" = "#FEAF16",
#                 "L3~Sp09D02" = "#00BCF9",
#                 "L4.inhib~Sp09D09" = "#C82100",
#                 "L5~Sp09D03" = "#16FF32",
#                 "L6~Sp09D04" = "#116A52",
#                 "WM.uf~Sp09D07" = "#E4E1E3",
#                 "WM~Sp09D06" = "#581009")

#### Spot plots for representative sections ####

single_vis_clus <- vis_clus(
    spe = spe,
    point_size = 1.5,
    colors = SpD_colors,
    sampleid = "Br5517",
    clustervar = "SpD",
    spatial = FALSE
) + 
    guides(fill = guide_legend(override.aes = list(size = 5)))

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

ggsave(cluster_plots[[1]][[1]], filename = ggsave(here(plot_dir, "test09.png")), width = 18)

cluster_grid <- Reduce("/", cluster_plots)
ggsave(cluster_grid, filename = here(plot_dir, "vis_SpD_rep_sections.pdf"), width = 18, height = 9)


#### Plot reduced dims ####
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(spe,
                          prefix = "ERC_spe",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "SpD",
                          color_pal = SpD_colors))


#### Spatial Registration vs. DLPFC ####
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


## Annotate modeling results
erc_modeling <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm", "modeling_results-BayesSpace_SVGm_k09.rds"))

erc_spd_anno <- readxl::read_excel(here("processed-data","05_spe_correct_cluster", "10_spatial_registration_DLPFC", "ERC_SpD_spatial_registration_Annotations.xlsx")) |>
    mutate(Annotation = fct_reorder(Annotation, order)) |>
    mutate(SpD = fct_reorder(paste0(Annotation, "~", cluster), order)) |>
    select(cluster, Annotation, SpD)

erc_modeling <- map(erc_modeling, function(mod){
    erc_colnames <- colnames(mod)
    pwalk(erc_spd_anno, function(...) erc_colnames <<- gsub(..1, ..3, erc_colnames))
    colnames(mod) <- erc_colnames
    return(mod)
})

map(erc_modeling, colnames)

saveRDS(erc_modeling, file = here(data_dir, "modeling_results-BayesSpace_SVGm_k09_annotated.rds"))

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
    confidence_threshold = 0.6,
    cutoff_merge_ratio = 0.25
))

## save data
save(cor_layer, anno, file = here(data_dir, "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"))

## refrence colors
layer_colors <- list(HumanPilot = spatialLIBD::libd_layer_colors,
                     spatialDLPFC = NULL) #TODO add spatial domain colors)
                     
map(names(cor_layer), function(ref){
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s.pdf", ref)))
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
    pdf(here(plot_dir, sprintf("layer_stat_cor_%s_cluster.pdf", ref)))
    print(layer_stat_cor_plot(
        cor_stats_layer = cor_layer[[ref]],
        reference_colors = layer_colors[[ref]],
        annotation = anno[[ref]],
        query_colors = SpD_colors
    ))
    dev.off()
})


# slurmjobs::job_single('22_SpD_clean_plots', create_shell = TRUE, memory = '25G', command = "Rscript 22_SpD_clean_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
