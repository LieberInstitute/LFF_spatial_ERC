## Adapted from 36_sn_subcluster_dotplot.R
## Plot NE (noradrenergic) receptor gene expression across cell types

library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("HDF5Array")
library("scDotPlot")
library("DeconvoBuddies")
library("tidyverse")

data_dir <- here("processed-data", "26_NE_exploration", "01_NEreceptor_expression")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- here("plots", "26_NE_exploration", "01_NEreceptor_expression")
if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))

rownames(sce) <- rowData(sce)$gene_name

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
load(here("processed-data","00_project_prep","Oligo_OPC_colors.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors

#### NE receptor genes ####

NE_receptor_genes <- c("ADRA1A","ADRA1B","ADRA1D","ADRA2A","ADRA2B","ADRA2C","ADRB1","ADRB2","ADRB3")

all(NE_receptor_genes %in% rownames(sce))
NE_receptor_genes[!NE_receptor_genes %in% rownames(sce)]

#### Plot dotplots + violin plots for a given cell subset ####

plot_NEreceptor_expression <- function(sce_subset, subset_label, color_pal) {

    ## scDotPlots - scaled and raw logcounts versions
    walk(c(TRUE, FALSE), function(scale_expr) {
        dotplot_NE <- sce_subset |>
            scDotPlot(features = NE_receptor_genes,
                      group = "cell_type_anno",
                      groupAnno = "cell_type_anno",
                      scale = scale_expr,
                      annoColors = list("cell_type_anno" = color_pal),
                      clusterRows = FALSE,
                      groupLegends = FALSE)

        suffix <- if (scale_expr) "" else "_logCounts"

        ggsave(dotplot_NE, filename = here(plot_dir, sprintf("sn_NEreceptor_dotplot_%s%s.pdf", subset_label, suffix)))
        ggsave(dotplot_NE, filename = here(plot_dir, sprintf("sn_NEreceptor_dotplot_%s%s.png", subset_label, suffix)))
    })

    ## violin plot
    plot_marker_express_List(
        sce = sce_subset,
        gene_list = list(NE_receptors = NE_receptor_genes),
        pdf_fn = here(plot_dir, sprintf("sn_NEreceptor_violin_%s.pdf", subset_label)),
        cellType_col = "cell_type_anno",
        gene_name_col = "gene_name",
        color_pal = color_pal
    )
}

#### All cell types ####
plot_NEreceptor_expression(sce, "all_cell_types", cell_type_colors$anno)

#### Astrocyte subtypes ####
plot_NEreceptor_expression(sce[, sce$cell_type_broad == "Astro"], "Astro_subtypes", cell_type_colors$anno)

#### Oligo + OPC ####
plot_NEreceptor_expression(sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC"], "OligoOPC", Oligo_OPC_colors)

#### Oligo subtypes only ####
plot_NEreceptor_expression(sce[, sce$cell_type_broad == "Oligo"], "Oligo_subtypes", Oligo_OPC_colors)

# slurmjobs::job_single('01_NEreceptor_expression', create_shell = TRUE, memory = '10G', command = "Rscript 01_NEreceptor_expression.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
