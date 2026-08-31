## Louise Huuki-Myers, May 2025
## Add sub-cluster data to sce object, explore and annotate clusters

library("SingleCellExperiment")
library("jaffelab")
library("tidyverse")
library("HDF5Array")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
library("readxl")

## source reduced dims function
source(here("code", "utils", "my_plot_reduced_dim.R"))

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "28_subcluster_update_sce")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "28_subcluster_update_sce")
if(!dir.exists(data_dir)) dir.create(data_dir)

## Load HD5F sce
message(Sys.time(), "- load Harmony corrected sce")
sce <- loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))
sce

colnames(colData(sce))

#### Update Ancestry data ####
## updated FLARE data from Feb 2025
load(here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.Rdata"), verbose = TRUE)

samples_ancestry <- samples_ancestry |> column_to_rownames("BrNum")

sce$Anc_Afr <- samples_ancestry[sce$sample_id, "Afr"]
sce$Anc_Eur <- samples_ancestry[sce$sample_id, "Eur"]

#### Add glia sub-clustering data ####

glia_cell_types <- c("Astro", "Endo", "Micro", "Oligo", "OPC")
names(glia_cell_types) <- glia_cell_types

sce$glia_subcluster <- NA

glia_subcluster_info <- map_dfr(glia_cell_types, function(ct){
    ## load clusters
    cluster_fn <- here("processed-data", "04_snRNA-seq", "25_sn_cluster_subtype", ct, sprintf("walktrap_snn_k10_subclusters_%s.Rdata", ct))
    load(cluster_fn)
    
    ## load cluster info
    ct_cluster_info <- read_csv(here("processed-data", "04_snRNA-seq", "26_sn_subtype_check", ct , sprintf("ERCsn_subtype_clustering_info-%s_k10.csv", ct)),
                                show_col_types = FALSE) |>
        mutate(cell_type_broad = ct)
    
    anno_table <- ct_cluster_info |> select(k10, cell_type_k10) |> column_to_rownames("k10")
    
    sce[,sce$cell_type_broad == ct]$glia_subcluster <<- anno_table[paste0("k", clusters),]
    
    return(ct_cluster_info)
})

table(sce$glia_subcluster)
table(sce$glia_subcluster, sce$cell_type_class)

#### Update cell_type_anno ####

sce$cell_type_anno_r1 <- sce$cell_type_anno
sce$cell_type_broad_r1 <- sce$cell_type_broad

## unfactor
sce$cell_type_anno <- as.character(sce$cell_type_anno)
sce$cell_type_broad <- as.character(sce$cell_type_broad)

## read in sub-cluster annotations
subcluster_anno <- read_excel(here("processed-data", "04_snRNA-seq", "27_sn_gila_check","ERCsn_glia_subcluster_info_anno.xlsx"))

## update neuron annotations in review (Aug 2026)
subcluster_neuron_anno <- read_excel(here("processed-data", "04_snRNA-seq", "27.5_sn_neuron_check","ERC_neuron_subcluster_details_annotated.xlsx"))

table(subcluster_anno$cell_type_anno, subcluster_anno$cell_type_broad)

anno_table <- subcluster_anno |> 
    select(cell_type_k10, cell_type_broad, cell_type_anno, passALL_metricQC) |>
    column_to_rownames("cell_type_k10")

anno_table <- anno_table[sce[,sce$cell_type_class == "glia"]$glia_subcluster,]
dim(anno_table)
head(anno_table)

nrow(anno_table) == ncol(sce[,sce$cell_type_class == "glia"])

## Update colData - cell colors
sce[,sce$cell_type_class == "glia"]$cell_type_broad <- anno_table$cell_type_broad
sce[,sce$cell_type_class == "glia"]$cell_type_anno <- anno_table$cell_type_anno
## Update passALL_metricQC
sce[,sce$cell_type_class == "glia"]$passALL_metricQC <- anno_table$passALL_metricQC

## apply the updated neuron annotations (in review, Aug 2026), same
## structure as the glia update above, but keyed on cell_type_anno itself
## (matching sce$cell_type_anno_r1) rather than a raw cell_type_k10 subcluster
## id - so cell_type_anno names are unchanged for neurons; this table updates
## cell_type_broad and passALL_metricQC for the existing names
anno_table_neuron <- subcluster_neuron_anno |>
    select(cell_type_anno, cell_type_broad, passALL_metricQC) |>
    column_to_rownames("cell_type_anno")

anno_table_neuron <- anno_table_neuron[sce[,sce$cell_type_class == "neuron"]$cell_type_anno_r1,]
dim(anno_table_neuron)
head(anno_table_neuron)

nrow(anno_table_neuron) == ncol(sce[,sce$cell_type_class == "neuron"])
stopifnot(!anyNA(anno_table_neuron$cell_type_broad))

sce[,sce$cell_type_class == "neuron"]$cell_type_broad <- anno_table_neuron$cell_type_broad
sce[,sce$cell_type_class == "neuron"]$passALL_metricQC <- anno_table_neuron$passALL_metricQC

table(sce$cell_type_broad, sce$cell_type_broad_r1)
table(sce$cell_type_anno, sce$cell_type_anno_r1)

table(sce$passALL_metricQC)
# FALSE   TRUE 
# 2717 122966

## Fix Pax6 broad cell type 
sce$cell_type_broad[sce$cell_type_anno == "Inhib.Pax6"] <- "Inhib"

table(sce$cell_type_anno, sce$cell_type_broad)

## update factors
cell_type_levels_broad <- colData(sce) |>
    as.data.frame() |> 
    count(cell_type_class, cell_type_broad) |>
    arrange(cell_type_class, cell_type_broad) |>
    pull(cell_type_broad)

sce$cell_type_broad <- factor(sce$cell_type_broad, levels =  cell_type_levels_broad)

cell_type_level_tb <- colData(sce) |>
    as.data.frame() |> 
    count(cell_type_class, cell_type_broad, cell_type_anno) |>
    filter(!grepl("QC|ambig", cell_type_anno)) |>
    arrange(cell_type_broad) 

# cat(cell_type_level_tb$cell_type_anno, sep = "\n")

#### Define colors for fine cell types ####

# load colors
load(here("processed-data","00_project_prep","cell_type_colors.V2.Rdata"), verbose = TRUE)


#### summary by cluster ####
message(Sys.time(), " - Summarize clusters")

pd <- as.data.frame(colData(sce))

cluster_info <- pd |>
    group_by(cell_type_class, cell_type_broad, cell_type_anno) |>
    summarize(n = n(),
              prop = n/ncol(sce),
              median_sum = median(sum),
              median_detected = median(detected),
              median_Mito_percent = median(subsets_Mito_percent),
              median_scDblFinder.score = median(scDblFinder.score),
              passALL_metricQC = all(passALL_metricQC)
    )

cluster_info |> arrange(median_sum)
cluster_info |> filter(!passALL_metricQC)

cluster_info |> filter(passALL_metricQC)

write.csv(cluster_info, file = here(data_dir, "ERC_sn_subcluster_info.csv"), row.names = FALSE)

#### Plot Quality Metrics ####

cell_class_cutoffs <- read.csv(here("processed-data", "04_snRNA-seq", "09_cluster_QC","ERC_sn_cell_class_cutoffs.csv"))

qc_violin_plot_all <- pd |> 
    select(cell_type_class, cell_type_broad, cell_type_anno, passALL_metricQC, sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(cell_type_anno, cell_type_broad, cell_type_class,passALL_metricQC), names_to = "metric") |>
    mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))) |>
    ggplot() +
    geom_violin(aes(x = cell_type_anno, y = value, fill = cell_type_anno, colour = passALL_metricQC), 
                scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
    scale_fill_manual(values = cell_type_colors$anno, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    geom_hline(data = cell_class_cutoffs, aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    geom_label(data = cell_class_cutoffs, aes(y = cutoff, label = cutoff_anno), x = 3, color = "blue", vjust = -.5) +
    facet_grid(metric~cell_type_class, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Quality Metric Value", x = "preliminary cell type cluster, snn k=10")

ggsave(qc_violin_plot_all, filename = here(plot_dir, "ERC_sn_subcluster_QCmetricViolin_ALL.png"), width = 12, height = 12)

## glia only
cell_class_cutoffs_g <- cell_class_cutoffs |> 
    filter(cell_type_class == "glia") |> 
    mutate(metric = ifelse(metric == "sum", "log10_sum", metric),
           cutoff = ifelse(metric == "log10_sum", log10(cutoff), cutoff))

qc_violin_plot_all_g <- pd |> 
    filter(cell_type_class == "glia") |>
    mutate(log10_sum = log10(sum)) |>
    select(cell_type_class, cell_type_broad, cell_type_anno, passALL_metricQC, log10_sum, detected, subsets_Mito_percent, scDblFinder.score)  |>
    pivot_longer(!c(cell_type_anno, cell_type_broad, cell_type_class,passALL_metricQC), names_to = "metric") |>
    # mutate(metric = factor(metric, levels = c("sum", "detected", "subsets_Mito_percent", "scDblFinder.score"))) |>
    ggplot() +
    geom_violin(aes(x = cell_type_anno, y = value, fill = cell_type_anno, colour = passALL_metricQC), 
                scale = "width", draw_quantiles = c(.25, 0.5, .75)) +
    scale_fill_manual(values = cell_type_colors$anno, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    theme_bw() +
    geom_hline(data = cell_class_cutoffs_g, aes(yintercept = cutoff), color = "blue", linetype = "dashed") +
    facet_grid(metric~cell_type_broad, scales = "free", space="free_x") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Quality Metric Value", x = "Glia Subcluster Cell Types")

ggsave(qc_violin_plot_all_g, filename = here(plot_dir, "ERC_sn_subcluster_QCmetricViolin_ALL_glia.png"), width = 12, height = 12)

#### n cell barplot ####
n_cell_barplot <- cluster_info |>
    ggplot(aes(x = cell_type_anno, y = n, fill = cell_type_anno, color = passALL_metricQC)) +
    geom_col() +
    geom_text(aes(label = n), size = 2, vjust = -0.5) +
    scale_fill_manual(values = cell_type_colors$anno, guide = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red")) +
    # facet_wrap(~cell_type_broad) +
    facet_grid(.~cell_type_broad, scales = "free_x", space="free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(n_cell_barplot, filename = here(plot_dir, "ERC_sn_subcluster_barplot_ct_n_qc.png"), width = 10)

#### Drop problamatic cluster ####
sce <- sce[,sce$passALL_metricQC]
sce <- sce[,!grepl("ambig", sce$cell_type_anno)]
dim(sce)

sce$cell_type_fine <- droplevels(sce$cell_type_fine)
sce$cell_type_anno <- factor(sce$cell_type_anno, levels = cell_type_level_tb$cell_type_anno)

levels(sce$cell_type_anno)
pd <- as.data.frame(colData(sce))

#### plot on UMAP + TSNE ####
message(Sys.time(), " - Plot UMAP + TSNE")

## plot annotated clusters
walk(c("UMAP", "TSNE"),
     ~my_plot_reduced_dim(sce,
                          prefix = "ERC_sn_subcluster",
                          var_type = "cat",
                          dimred = .x,
                          my_var = "cell_type_anno",
                          color_pal = cell_type_colors$anno,
                          suffix = "sctype"))

#### plot marker genes ####
message(Sys.time(), " - Plot marker genes")

## read in marker genes from lit
lit_markers <- read_csv(here("processed-data","04_snRNA-seq", "00_lit_marker_genes", "lit_marker_summary.csv")) |>
    mutate(in_data = gene_name %in% rowData(sce)$Symbol,
           cell_type_broad = factor(cell_type_broad, levels = levels(sce$cell_type_broad))) |>
    arrange(cell_type_broad)

## missing from our data
lit_markers |> filter(!in_data)
# gene_name cell_type          n_studies studies                  cell_type_broad in_data
# <chr>     <chr>                  <dbl> <chr>                    <chr>           <lgl>  
# 1 CEMP      Fibroblast                 1 Davila-Velderrain et al. Endo            FALSE  
# 2 GR1A1     Glutamate receptor         1 Grubman et al.           Excit           FALSE  
# 3 PDCH15    OPC                        1 Grubman et al.           OPC             FALSE 

lit_markers <- lit_markers |> filter(in_data)

lit_markers_list <- map(splitit(lit_markers$cell_type), ~lit_markers$gene_name[.x])

lit_marker_order <- lit_markers |> count(cell_type_broad, cell_type) |> pull(cell_type)
lit_markers_list <- lit_markers_list[lit_marker_order]

plot_marker_express_List(sce, 
                         lit_markers_list, 
                         pdf_fn = here(plot_dir, "ERC_sn_sctype_lit_markers.pdf"),
                         cellType_col = "cell_type_anno",
                         gene_name_col = "gene_name",
                         color_pal = cell_type_colors$anno,
)


####  Cell Type Proportions ####
cell_class_proportions <- pd |>
    group_by(sample_id, APOE, Sex, Ancestry, cell_type_class) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n))

sample_neuron_order <- cell_class_proportions |> filter(cell_type_class == "neuron") |> arrange(prop) |> pull(sample_id)

cell_type_proportions <- pd |>
    group_by(sample_id, APOE, APOE_carrier, Sex, Age, Ancestry, Anc_Afr, exp_round, seq_round, cell_type_anno) |>
    summarize(n = n()) |> 
    group_by(sample_id, APOE, Sex, Ancestry) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = sample_neuron_order))

## save proportion data
save(cell_type_proportions, file = here(data_dir, "cell_type_proportions.Rdata"))
write.csv(cell_type_proportions, file = here(data_dir, "cell_type_proportions.csv"))

## plot n per sample
cell_type_sample_n_bar <- cell_type_proportions |>
    ggplot(aes(x = sample_id, y = n, fill = cell_type_anno)) +
    geom_col() +
    geom_text(aes(label = ifelse(n > 100, n, "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors$anno) +
    labs(y = "n Nucei")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    guides(fill=guide_legend(ncol=1))

ggsave(cell_type_sample_n_bar, filename = here(plot_dir, "ERC_sn_subtype_barplot_ct_sample_n.png"), width = 10)

## plot proprotions
cell_type_proportion_bar <- cell_type_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = cell_type_anno)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors$anno) +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    guides(fill=guide_legend(ncol=1))

ggsave(cell_type_proportion_bar, filename = here(plot_dir, "ERC_sn_subtype_barplot_ct_prop.png"), width = 10)

cell_type_proportion_bar_facet <- cell_type_proportions |>
    ggplot(aes(x = sample_id, y = prop, fill = cell_type_anno)) +
    geom_col() +
    geom_text(aes(label = ifelse(prop > .02, round(prop, 2), "")),
              position = position_stack(vjust = .5),
              size = 2) +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors$anno) +
    labs(y = "Cell Type Proportion")  +
    facet_grid(.~APOE, scales = "free_x", space = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    guides(fill=guide_legend(ncol=1)) 

ggsave(cell_type_proportion_bar_facet, filename = here(plot_dir, "ERC_sn_subtype_barplot_ct_prop_facet.png"), width = 10)

## evaluate proportion sample by cell 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
# we need sample_colors

sample_proportions <- pd |>
    group_by(sample_id, cell_type_anno) |>
    summarize(n = n()) |> 
    group_by(cell_type_anno) |>
    mutate(prop = n/sum(n),
           sample_id = factor(sample_id, levels = names(sample_colors))) 

## check for cell types from few samples 
sample_proportions |> arrange(-prop)

sample_qc <- sample_proportions |> 
    group_by(cell_type_anno) |>
    summarize(n_samples = length(unique(sample_id)),
              max_prop = max(prop)) |>
    mutate(pass_sampleQC = max_prop < 0.5 & n_samples > 5)

## lowest contibution is 19 samples 
sample_qc |> arrange(n_samples)
sample_qc |> arrange(-max_prop)

## n donor by cell type
ct_n_donor <- sample_qc |>
    ggplot(aes(y = reorder(cell_type_anno, n_samples), x = n_samples, fill = cell_type_anno)) +
    geom_col() +
    geom_text(aes(label = n_samples), hjust = 1) +
    scale_fill_manual(values = cell_type_colors$anno) +
    theme_bw() +
    labs(x = "n Donors", y = "Cell Type") +
    theme(legend.position = "None")

ggsave(ct_n_donor, filename = here(plot_dir, "ERC_sn_subtype_barplot_ct_n_donor.png"), width = 5, height = 6)


## plot sample by cell type proportions bar
sample_n_bar <- sample_proportions |>
    ggplot(aes(x = cell_type_anno, y = n, fill = sample_id)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values = sample_colors) +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sample_n_bar, filename = here(plot_dir, "ERC_sn_subtype_barplot_sample_n.png"), width = 10)

sample_proportion_bar <- sample_proportions |>
    ggplot(aes(x = cell_type_anno, y = prop, fill = sample_id)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values = sample_colors) +
    labs(y = "Cell Type Proportion")  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sample_proportion_bar, filename = here(plot_dir, "ERC_sn_subtype_barplot_sample_prop.png"), width = 10)

#### Add colors to metadata ####
cell_type_colors <- list(broad = cell_type_colors$broad[levels(sce$cell_type_broad)],
                         anno = cell_type_colors$anno[levels(sce$cell_type_anno)])

metadata(sce)$cell_type_colors <- cell_type_colors

#### Output annotations ####
message(Sys.time(), " - Save annotated sce")

# save(sce, file = here("processed-data", "spe_objects", "sce_ERC.Rdata"))
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce_objects", "sce_ERC_subcluster"), replace=TRUE)

# slurmjobs::job_single('28_subcluster_update_sce', create_shell = TRUE, memory = '10G', command = "Rscript 28_subcluster_update_sce.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

