## Louise Huuki-Myers, Aug 2026
## Gather neuron details (pre-print version) to help re-annotate some labels

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
# library("HDF5Array")
library("spatialLIBD")
library("tidyverse")
library("jaffelab")
library("DeconvoBuddies")
library("bluster")
library("ComplexHeatmap")

data_dir <- here("processed-data", "04_snRNA-seq", "27.5_sn_neuron_check")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# plot_dir <- here("plots", "04_snRNA-seq", "27.5_sn_neuron_check")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### general info ####

## first round clustering notes
anno_notes <- read_csv(here("processed-data", "04_snRNA-seq", "15_cluster_update_sce", "ERC_sn_cluster_info.csv")) |> select(snn_kOpt, cell_type_anno)

subcluster_info <- read_csv(here("processed-data", "04_snRNA-seq", "28_subcluster_update_sce", "ERC_sn_subcluster_info.csv"))|>
    filter(cell_type_class == "neuron") |>
    left_join(anno_notes)

#### registration annotations ####

## xenium density
cell_v_xSpD_prop <- read_csv(here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "cell_v_xSpD_prop_long.csv")) |>
    select(-`...1`)

cell_v_xSpD_prop |> filter(cell_type_anno == "Inhib.Lamp5_Lhx6")

neuron_v_xSpD_density <- cell_v_xSpD_prop |> 
    filter(grepl("Excit|Inhib", cell_type_anno), prop_xSpD > 0.19) |>
    group_by(cell_type_anno) |>
    summarise(xSpD = paste0(xSpD, " (", round(100 * prop_xSpD), "%)", collapse = ", "))
    

## regular registration 
registration_anno_summary <- read_csv(here("processed-data", "04_snRNA-seq", "30_sn_subcluster_spatial_registration_anno", "ERCsn_subcluster_spatial_registration_anno_summary.csv")) |> 
    dplyr::rename(cell_type_anno = cluster) |>
    select(cell_type_anno, spatialERC, snDLPFC_PEC, sestan_EC, spatialHPC) |>
    filter(grepl("Excit|Inhib", cell_type_anno)) |>
    left_join(neuron_v_xSpD_density)

registration_anno_summary

#### Marker genes ####

# marker_stats_MeanRatio - global
# load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio", "marker_stats_MeanRatio_cell_type_anno.Rdata"), verbose = TRUE)
# 
# marker_stats_MeanRatio |> 
#     filter(grepl("Excit|Inhib", cellType.target), !grepl("^ENSG", gene), MeanRatio > 1, MeanRatio.rank <=10)

## subtype specific 
top_enrich_genes <- map_dfr(c("Excit", "Inhib"), ~read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", .x, sprintf("subtype_enrichment_top10_%s.csv", .x)))) |>
    filter(!grepl("^ENSG", gene))

top_enrich_genes_list <- top_enrich_genes |>
    group_by(cell_type_anno = test) |>
    summarize(enrichment_genes = paste(gene, collapse = ", "))

top_MeanRatio_genes <- map_dfr(c("Excit", "Inhib"), ~read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", .x, sprintf("subtype_MeanRatio_top10_%s.csv", .x)))) |>
    filter(!grepl("^ENSG", gene), MeanRatio > 1)

top_MeanRatio_genes_list <- top_MeanRatio_genes |>
    group_by(cell_type_anno = cellType.target) |>
    summarize(MeanRatio_genes = paste(gene, collapse = ", "))

#### Combine data ####

neuron_subcluster_details <- subcluster_info |>
    left_join(registration_anno_summary) |>
    left_join(top_enrich_genes_list) |>
    left_join(top_MeanRatio_genes_list)

write_csv(neuron_subcluster_details, file = here(data_dir, "ERC_neuron_subcluster_details_preprint_label.csv"))

# slurmjobs::job_single(name = "27_sn_glia_check", memory = "50G", command = "Rscript 27_sn_glia_check.R", create_shell = TRUE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



