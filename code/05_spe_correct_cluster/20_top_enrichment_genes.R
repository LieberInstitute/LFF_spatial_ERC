## Dec 2024, Louise Huuki-Myers
## Extract and plot top enrichment model genes from selected cluster

library("spatialLIBD")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "05_spe_correct_cluster", "20_top_enrichment_genes")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "05_spe_correct_cluster", "20_top_enrichment_genes")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load SPE data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

## load psuedobulk and modeling
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","spe_pseudobulk-BayesSpace_SVGm_k09.rds"))
modeling_results <- readRDS(here("processed-data", "05_spe_correct_cluster", "08_model_pseudobulk", "BayesSpace_SVGm","modeling_results-BayesSpace_SVGm_k09.rds"))

## annotation
load(here("processed-data", "05_spe_correct_cluster", "10_spatial_registration_DLPFC", "spatial_registration_erc_v_DLPFC_cor_anno.Rdata"), verbose = TRUE)
# cor_anno$BayesSpace_SVGm_k09$layer_anno$HumanPilot

cluster_anno_k9 <- cor_anno$BayesSpace_SVGm_k09$layer_anno$HumanPilot |>
    mutate(HumanPilot_anno = paste0(layer_label_simple, "~", cluster)) |>
    select(test = cluster, HumanPilot_anno, layer_label = layer_label_simple)

#### Extract Top Layer Enrichment Genes ####
top_DEGs <- sig_genes_extract(n = 100,
                              modeling_results = modeling_results,
                              model_type = "enrichment",
                              sce_layer = spe_pb) |>
    left_join(cluster_anno_k9)

write_csv(top_DEGs, file = here(data_dir, "enrichment_modeling_top100-BayesSpace_SVGm_k09.csv"))

#### Top Enrichment Gene Heatmap ####
top10_DEGs <- top_DEGs |>
    filter(top <= 10) |>
    arrange(HumanPilot_anno)

## extract and scale logcounts
rownames(spe_pb) <- rowData(spe_pb)$gene_name
logcounts_z <- t(scale(t(logcounts(spe_pb)[top10_DEGs$gene,])))

dim(logcounts_z)
# [1]  90 279


## colors 
# L1~Sp09D05   L1~Sp09D08 L2/3~Sp09D01   L3~Sp09D02   L4~Sp09D09 L5/4~Sp09D03   L6~Sp09D04   WM~Sp09D06   WM~Sp09D07
k09_colors <- c(Sp09D05 = "#16FF32",
                Sp09D08 = "#90AD1C",
                Sp09D01 = "#5A5156",
                Sp09D02 = "#3283FE",
                Sp09D09 = "brown",
                Sp09D03 = "#FE00FA",
                Sp09D04 =  "#F6222E",
                Sp09D06 = "#FEAF16",
                Sp09D07 = "#1CFFCE")

layer_label_colors <- spatialLIBD::libd_layer_colors
names(layer_label_colors) <- gsub("ayer", "", names(layer_label_colors))
layer_label_colors <- c(layer_label_colors, "L2/3" = "#50DDAC", "L5/4" = "#BD8339")
layer_label_colors <- layer_label_colors[unique(top10_DEGs$layer_label)]
# L2/3        L3      L5/4        L6        L1        WM        L4 
# "#50DDAC" "#4DAF4A" "#BD8339" "#FF7F00" "#F0027F" "#1A1A1A" "#984EA3" 

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

## pb col annotations
pb_pd <- as.data.frame(colData(spe_pb)) 

pb_anno <- pb_pd |>
    select(#n_spots = ncells,
           # sample_id,
           Sex,
           APOE,
           Ancestry,
           k09 = BayesSpace_SVGm_k09)

pb_column_ha <- HeatmapAnnotation(
   df = pb_anno,
   col = list(k09 = k09_colors,
              Sex = sex_colors, 
              APOE = APOE_genotype_colors,
              Ancestry = ancestry_colors)
)


## gene row annotations
gene_anno <- top10_DEGs |>
    select(k09 = test,
           layer_label)

unique(top10_DEGs$layer_label)

gene_row_ha <- rowAnnotation(df = gene_anno,
                             col = list(k09 = k09_colors,
                                        layer_label = layer_label_colors))

## plot heatmap
pdf(here(plot_dir, "enrichment_top10_heatmap-BayesSpace_SVGm_k09.pdf"), height = 15, width = 10)
Heatmap(logcounts_z,
        right_annotation = gene_row_ha,
        top_annotation = pb_column_ha,
        cluster_rows = FALSE,
        # cluster_columns = FALSE,
        show_column_names = FALSE,
        row_split = gene_anno$k09
        )
dev.off()


#### Spot plots ####



# slurmjobs::job_single('20_top_enrichment_genes', create_shell = TRUE, memory = '25G', command = "Rcript 20_top_enrichment_genes.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
