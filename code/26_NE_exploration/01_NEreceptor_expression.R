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

plot_NEreceptor_expression <- function(sce_subset, subset_label, color_pal, cell_type_col = "cell_type_anno") {

    my_colors = list(color_pal)
    names(my_colors) <- cell_type_col
    
    ## scDotPlots - scaled and raw logcounts versions
    walk(c(FALSE, TRUE), function(scale_expr) {
        dotplot_NE <- sce_subset |>
                scDotPlot(features = NE_receptor_genes,
                          group = cell_type_col,
                          groupAnno = cell_type_col,
                          scale = scale_expr,
                          annoColors = my_colors,
                          clusterRows = FALSE,
                          clusterColumns = FALSE,
                          groupLegends = FALSE)
            
            suffix <- ifelse(scale_expr,  "" , "_logCounts")
            
            ggsave(dotplot_NE, filename = here(plot_dir, sprintf("sn_NEreceptor_dotplot_%s%s.pdf", subset_label, suffix)))
            ggsave(dotplot_NE, filename = here(plot_dir, sprintf("sn_NEreceptor_dotplot_%s%s.png", subset_label, suffix)))
            })

    gene_violin <- plot_gene_express(
        sce = sce_subset,
        genes = NE_receptor_genes,
        category = cell_type_col,
        color_pal = color_pal,
        free_y = TRUE
    )
    
    ggsave(gene_violin, filename = here(plot_dir, sprintf("sn_NEreceptor_violin_%s.pdf", subset_label)))
    ggsave(gene_violin, filename = here(plot_dir, sprintf("sn_NEreceptor_violin_%s.png", subset_label)))
    
    
    ADRA1A_violin <- plot_gene_express(
        sce = sce_subset,
        genes = "ADRA1A",
        category = cell_type_col,
        color_pal = color_pal,
        free_y = TRUE
    )
    
    ggsave(ADRA1A_violin, filename = here(plot_dir, sprintf("sn_ADRA1A_violin_%s.png", subset_label)), height = 4, width = 6)
}

#### All cell types ####
plot_NEreceptor_expression(sce, "all_cell_types", cell_type_colors$anno)

#### Astrocyte subtypes ####
plot_NEreceptor_expression(sce_subset <- sce[, sce$cell_type_broad == "Astro"], 
                           subset_label = "Astro_subtypes", 
                           color_pal <- cell_type_colors$anno)

#### Oligo + OPC ####
plot_NEreceptor_expression(sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC"], "OligoOPC", Oligo_OPC_colors)

#### Oligo subtypes only ####
plot_NEreceptor_expression(sce[, sce$cell_type_broad == "Oligo"], "Oligo_subtypes", Oligo_OPC_colors)

#### Oligo annotated grouos ####
sce_oligo <- sce[, sce$cell_type_broad == "Oligo" | sce$cell_type_broad == "OPC"]
sce_oligo$cell_type_anno <- droplevels(sce_oligo$cell_type_anno)

oligo_notes <- readxl::read_xlsx(here("processed-data", "19_other_Oligo", "00_check_Oligo_markers", "ERC_Oligo_notes.xlsx")) |>
    select(ERC_Oligo = `Oligo subtype`, Oligo_anno = Annotation_short) |>
    add_row(ERC_Oligo = paste0("OPC.", 1:5), Oligo_anno = "OPC") |>
    mutate(Oligo_anno = fct_inorder(Oligo_anno))

levels(oligo_notes$Oligo_anno)

sce_oligo$Oligo_anno <- oligo_notes$Oligo_anno[match(sce_oligo$cell_type_anno, oligo_notes$ERC_Oligo)]

table(sce_oligo$Oligo_anno, sce_oligo$cell_type_anno)

Oligo_anno_colors <- c(Oligo.M  = "#00C4B8",
                       Oligo.PM  = "#00767A", 
                       Oligo.L  = "#0072E5",  
                       Oligo.NF = "#7B61FF", 
                       Oligo.PV  = "#FF3D9A",  
                       OPC      = "#C49A00",
                       Oligo.other = "grey")

plot_NEreceptor_expression(sce_oligo, "Oligo_annotations", color_pal = Oligo_anno_colors, cell_type_col = "Oligo_anno")



# slurmjobs::job_single('01_NEreceptor_expression', create_shell = TRUE, memory = '10G', command = "Rscript 01_NEreceptor_expression.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
