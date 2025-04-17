## Louise Huuki-Myers, April 2025
## Examine covaraites in Visium dataset for DGE
## adapted from https://github.com/LieberInstitute/dlpfc_asd/blob/a250a1d7e20bd754c5f1186aa96ce0752d55e556/code/08_pseudoBulkDGE_s/02_covariate_analysis.R

library("spatialLIBD")
library("SingleCellExperiment")
# library("scran")
# library("BayesSpace")
library("tidyverse")
library("ggpubr")
library("ggrepel")
library("here")
library("sessioninfo")


#### Set up dirs ####
plot_dir <- here("plots", "09_pseudoBulkDGE_Visium", "02_covariate_analysis_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
data <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))

## load colors
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

# Extract the relevant columns from the data
plot_data <- as.data.frame(colData(data)) |>
    select(sample_id, SpD, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, Visium_slide, round, nspots = ncells, pseudo_sum_umi, pseudo_expr_chrM_ratio) |>
    mutate(APOE = factor(APOE))
           

## TODO add pseudobulk quality metrics sum_umi, expr_chrM_ratio, pmi

# Generate n cell boxplots for each cell_type

pdf(here(plot_dir, "ERC_Visium_nspots_vs_APOE.pdf"), width = 8, height = 6)
for (d in levels(plot_data$SpD)) {
    
    cell_type_data <- plot_data |> filter(SpD == d)
    y_max_ncells <- max(cell_type_data$nspots)*1.1
    
    plot <- ggboxplot(
        cell_type_data, 
        x = "APOE_carrier", 
        y = "nspots", 
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
        ggtitle(paste("SpD:", d)) +
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