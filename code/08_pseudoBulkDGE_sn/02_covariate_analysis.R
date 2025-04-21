## Louise Huuki-Myers, April 2025
## Examine covaraites in dataset for DGE
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

sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data","sce_pseudo_DGE.RDS"))

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

colnames(colData(sce_pb))

# Extract the relevant columns from the data
plot_data <- as.data.frame(colData(sce_pb)) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE,
           cell_type_anno, ncells, exp_round, seq_round,
           pseudo_sum_umi, pseudo_expr_chrM, pseudo_expr_chrM_ratio) 

levels(plot_data$APOE)

# Generate n cell boxplots for each cell_type

cell_type_count <- plot_data |> group_by(APOE_carrier, cell_type_anno) |> count()

pb_cell_type_bar_APOE_carrier <- plot_data |>
    count(APOE_carrier, cell_type_anno) |> 
    ggplot(aes(x = cell_type_anno, y=n, fill = APOE_carrier)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_carrier_colors)

ggsave(pb_cell_type_bar_APOE_carrier, filename = here(plot_dir, 'pb_cell_type_bar_APOE_carrier.png'))

pb_cell_type_bar_APOE <- plot_data |>
    count(APOE, cell_type_anno) |> 
    ggplot(aes(x = cell_type_anno, y=n, fill = APOE)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_genotype_colors)

ggsave(pb_cell_type_bar_APOE, filename = here(plot_dir, 'pb_cell_type_bar_APOE.png'))


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

#### variance parition on Sample level data ####
library("variancePartition")
# load sample level data

sce_pb_sample <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk", "sce_pseudobulk_only-sample_id.rds"))
dim(sce_pb_sample)
table(sce_pb_sample$APOE)
# E2/E2 E2/E3 E3/E4 E4/E4 
# 6     8    10     7 
table(sce_pb_sample$APOE_carrier)
# E2+ E4+ 
# 14  17 

pd_sample <- as.data.frame(colData(sce_pb_sample))

# Assess correlation between all pairs of variables
colnames(pd_sample)

# form <- ~ APOE + sample_id + Sex + Age + Anc_Afr + exp_round + seq_round + ncells
form <- ~ APOE + sample_id + exp_round
# form <- ~ APOE + Sex + Age + Anc_Afr + exp_round + seq_round + ncells

C <- canCorPairs(form, pd_sample)
# ~ APOE + sample_id + Sex + Age + Anc_Afr + exp_round + seq_round + ncells
# APOE and sample_id
# sample_id and Sex
# sample_id and Age
# sample_id and Anc_Afr
# sample_id and exp_round
# sample_id and seq_round
# sample_id and ncells

pdf(here(plot_dir, "sn_variable_cor_matrix.pdf"))
plotCorrMatrix(C)
dev.off()

## test variance partition 
form <- ~ (1 | APOE_carrier) + (1 | APOE) + (1 | Sex) + Anc_Afr + Age + (1 | exp_round) + (1 | seq_round)
varPart <- fitExtractVarPartModel(logcounts(sce_pb_sample), form, pd_sample)
vp <- sortCols(varPart)

## violin plot 
sn_vp_violin <- plotVarPart(vp)
ggsave(sn_vp_violin, filename = here(plot_dir, "sn_vp_violin.png"))

# slurmjobs::job_single('01_create_pseudobulk_data', create_shell = TRUE, memory = '25G', command = "Rscript 01_create_pseudobulk_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()