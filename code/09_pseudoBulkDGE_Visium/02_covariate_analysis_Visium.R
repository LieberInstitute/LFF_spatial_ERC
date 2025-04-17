## Louise Huuki-Myers, April 2025
## Examine covaraites in Visium dataset for DGE
## adapted from https://github.com/LieberInstitute/dlpfc_asd/blob/a250a1d7e20bd754c5f1186aa96ce0752d55e556/code/08_pseudoBulkDGE_s/02_covariate_analysis.R

library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
# library("BayesSpace")
library("tidyverse")
library("ggpubr")
library("ggrepel")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data")
#if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Set up dirs ####
plot_dir <- here("plots", "08_pseudoBulkDGE_sn", "02_covariate_analysis")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
# message(Sys.time(), " - Load HDF5 sce")
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

data <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk","sce_pseudobulk-cell_type_anno.rds"))

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

# Extract the relevant columns from the data
plot_data <- as.data.frame(colData(data)) |>
    select(sample_id, cell_type_anno, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, exp_round, seq_round, ncells) |>
    mutate(APOE = factor(gsub("^(E[2,3,4])(E[2,3,4])","\\1/\\2", APOE)))
           
levels(plot_data$APOE)

## TODO add pseudobulk quality metrics sum_umi, expr_chrM_ratio, pmi

# Generate n cell boxplots for each cell_type

pdf(here(plot_dir, "ERC_sn_ncells_vs_APOE.pdf"), width = 8, height = 6)
for (cell_type in levels(plot_data$cell_type_anno)) {
    
    cell_type_data <- plot_data |> filter(cell_type_anno == cell_type)
    y_max_ncells <- max(cell_type_data$ncells)*1.1
    
    plot <- ggboxplot(
        cell_type_data, 
        x = "APOE_carrier", 
        y = "ncells", 
        color = "APOE_carrier", 
        palette = APOE_carrier_colors,
        add = "jitter", 
        shape = 19, 
        xlab = "APOE Genotype", 
        ylab = "Number of Nuclei"
    ) + 
        geom_text_repel(
            aes(label = sample_id),
            position = position_jitter(width = 0.1, height = 0.2),
            hjust = 0.5, 
            vjust = 1,
            size = 3
        ) + 
        ggtitle(paste("Cell Type:", cell_type)) +
        theme_bw() + 
        theme(legend.position = "none") +
        stat_compare_means(
            aes(group = APOE_carrier, label = paste0("p = ", after_stat(p.format))),
                           label.x = 0.5,
                           label.y = y_max_ncells,
                           method = "t.test")
    print(plot)
}
dev.off()

## bar plot of n samples by cell type that pass pseudobulk cutoff 

# slurmjobs::job_single('01_create_pseudobulk_data', create_shell = TRUE, memory = '25G', command = "Rscript 01_create_pseudobulk_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()