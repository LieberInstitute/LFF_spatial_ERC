## Louise Huuki-Myers, Aug 2025
## plot cell type marker genes in Spatial data

#### set up ####

library("HDF5Array")
library("spatialLIBD")
library("sessioninfo")
library("here")
library("tidyverse")
library("patchwork")
library("escheR")

plot_dir <- here("plots", "15_spot_deconvolution", "04_ct_marker_gene_spatial_expression")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "15_spot_deconvolution", "04_ct_marker_gene_spatial_expression")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

#### load spatial data ####
message(Sys.time(), " - Load Spatial data")
## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
rownames(spe) <- rowData(spe)$gene_name

## load deconvolution data

spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("processed-data", "spe_objects", "spe_ERC_annotated"))
rownames(spe) <- rowData(spe)$gene_name


## load marker data

celltype = "Oligo"

load(here("processed-data", "04_snRNA-seq", "34_sn_subcluster_MeanRatio", "marker_stats_MeanRatio_cell_type_anno.Rdata"),
     verbose = TRUE)

MR_global_top25 <- marker_stats_MeanRatio |> 
    filter(grepl(celltype, cellType.target),
           MeanRatio > 1, 
           MeanRatio.rank <= 25) |>
    mutate(cellType.target = droplevels(cellType.target))

MR_global_top25 |> count(cellType.target)

MR_global_list <- map(rafalib::splitit(MR_global_top25$cellType.target), ~MR_global_top25[.x,] |> pull(gene))
MR_global_list <- map(MR_global_list, ~.x[.x %in% rownames(spe)])

## subtype markers 
list.files(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", celltype))

ct_enrich_top10 <- read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", celltype,
                        sprintf("subtype_enrichment_top10_%s.csv", celltype)), row.names = 1)

enrich_gene_list <- map(rafalib::splitit(ct_enrich_top10$test), ~ct_enrich_top10[.x, "gene"])

ct_MR_top10 <- read.csv(here("processed-data", "04_snRNA-seq", "35_sn_subcluster_marker_modeling", celltype,
                        sprintf("subtype_MeanRatio_top10_%s.csv", celltype)))

MR_gene_list <- map(rafalib::splitit(ct_MR_top10$cellType.target), ~ct_MR_top10[.x, "gene"])
MR_gene_list <- map(MR_gene_list, ~.x[.x %in% rownames(spe)])

#### Vis genes ####
source(here("code", "15_spot_deconvolution", "vis_rep_sections.R"))

test_vis_gene <- vis_gene(spe,
                          geneid = MR_global_list$Oligo.1,
                          assayname = "logcounts",
                          multi_gene_method = "z_score",
                          cont_colors = viridisLite::rocket(10, direction = -1))

ggsave(test_vis_gene, filename = here(plot_dir, "vis_gene_test.png"))

vis_rep_test <- vis_rep_sections(spe, geneid = MR_gene_list$Oligo.1)
ggsave(vis_rep_test, filename = here(plot_dir, "vis_rep_test.png"), width = 18, height = 9)


#### escheR plots ####

spe$counts_MOBP <- counts(spe)[which(rowData(spe)$gene_name=="MOBP"),]

p <- make_escheR(spe[, spe$sample_id %in% c("Br5517")]) |>
    add_ground(var = "SpD") |>
    add_fill(var = "counts_MOBP") + 
    scale_fill_gradient(low = "white", high = "black") +
    scale_color_manual(values = SpD_colors)

ggsave(p, filename = here(plot_dir, "escher_test.png"))

