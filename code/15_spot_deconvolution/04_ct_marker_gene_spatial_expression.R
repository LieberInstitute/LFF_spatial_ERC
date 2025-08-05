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
library("DeconvoBuddies")

plot_dir <- here("plots", "15_spot_deconvolution", "04_ct_marker_gene_spatial_expression")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "15_spot_deconvolution", "04_ct_marker_gene_spatial_expression")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

celltype = "Oligo"

#### load spatial data ####
message(Sys.time(), " - Load Spatial data")
## Load HD5F spe
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC_annotated"))
rownames(spe) <- rowData(spe)$gene_name

#### check RCTD data ####

## load deconvolution data
spe_rctd <- readRDS(here("processed-data", "15_spot_deconvolution", "03_results_RCTD","cell_type_broad", "spe_RCTD-cell_type_broad.rds"))

## add to rowData
summary(assay(spe_rctd, "weights")[celltype,])
table(assay(spe_rctd, "weights")[celltype,] > 0)

colData(spe)$RCTD_weight <- NA
spe[,colnames(spe_rctd)]$RCTD_weight <- assay(spe_rctd, "weights")[celltype,]

spe$RCTD_ct <- spe$RCTD >0
table(spe$SpD, spe$RCTD_ct)

spe$RCDT_SpD <- ifelse(spe$RCTD_ct, paste0(spe$SpD, "_", celltype), paste0(spe$SpD, "_no", celltype))
table(spe$RCDT_SpD)

RCDT_SpD_colors <- SpD_colors
names(RCDT_SpD_colors) <- paste0(names(RCDT_SpD_colors), "_", celltype)

RCDTno_SpD_colors <- SpD_colors
names(RCDTno_SpD_colors) <- paste0(names(RCDTno_SpD_colors), "_no", celltype)

RCDT_SpD_colors <- sort(c(RCDT_SpD_colors, RCDTno_SpD_colors))

#### load marker data ####
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

## plot marker data
plot_marker_express_List(
    spe,
    gene_list = enrich_gene_list,
    pdf_fn = here(plot_dir, sprintf("SpD_rctd_enrich-%s.pdf", celltype)),
    cellType_col = "RCDT_SpD",
    color_pal = RCDT_SpD_colors,
    plot_points = FALSE
)


#### Vis genes ####
source(here("code", "15_spot_deconvolution", "vis_rep_sections.R"))

# vis_rep_test <- vis_rep_sections(spe, geneid = MR_gene_list$Oligo.1)
# ggsave(vis_rep_test, filename = here(plot_dir, "vis_rep_test.png"), width = 18, height = 9)

pdf(here(plot_dir, sprintf("vis_rep_enrich_markers-%s.pdf", celltype)), width = 18, height = 9)
map(enrich_gene_list, ~vis_rep_sections(spe, geneid = .x))
dev.off()


#### escheR plots ####

spe$counts_MOBP <- counts(spe)[which(rowData(spe)$gene_name=="MOBP"),]

p <- make_escheR(spe[, spe$sample_id %in% c("Br5517")]) |>
    add_ground(var = "SpD") |>
    add_fill(var = "counts_MOBP") + 
    scale_fill_gradient(low = "white", high = "black") +
    scale_color_manual(values = SpD_colors)

ggsave(p, filename = here(plot_dir, "escher_test.png"))







