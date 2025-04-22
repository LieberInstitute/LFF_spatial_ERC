## Louise Huuki-Myers, April 2025
## Examine covairates in sn dataset for DGE
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
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "02_covariate_analysis")
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
pd <- as.data.frame(colData(sce_pb)) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE,
           cell_type_anno, ncells, exp_round, seq_round,
           pseudo_sum_umi, pseudo_expr_chrM, pseudo_expr_chrM_ratio) 


#### Generate n sample barplots for each cell_type ####

cell_type_count <- pd |> group_by(APOE_carrier, cell_type_anno) |> count()

pb_cell_type_bar_APOE_carrier <- pd |>
    count(APOE_carrier, cell_type_anno) |> 
    ggplot(aes(x = cell_type_anno, y=n, fill = APOE_carrier)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_carrier_colors)

ggsave(pb_cell_type_bar_APOE_carrier, filename = here(plot_dir, 'pb_cell_type_bar_APOE_carrier.png'))

pb_cell_type_bar_APOE <- pd |>
    count(APOE, cell_type_anno) |> 
    ggplot(aes(x = cell_type_anno, y=n, fill = APOE)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_genotype_colors)

ggsave(pb_cell_type_bar_APOE, filename = here(plot_dir, 'pb_cell_type_bar_APOE.png'))

#### t-test variables ####

## ncells 

y_position <- pd |>
    group_by(cell_type_anno) |>
    summarise(y.position = max(ncells) + .05*max(ncells))

n_cell_t_test <- pd |>
    filter(cell_type_anno != "Excit.L5_6_NP") |>
    mutate(cell_type_anno = droplevels(cell_type_anno)) |>
    do(compare_means(ncells ~ APOE_carrier, data = ., method = "t.test", p.adjust.method = "fdr", group.by = "cell_type_anno")) |>
    ungroup() |>
    mutate(p.signif.fdr = case_when(p.adj < 0.005 ~ "***",
                                     p.adj < 0.01 ~"**",
                                    p.adj < 0.05 ~"*",
                                     TRUE~""),
           fdr_anno = sprintf("FDR=%.3f%s", p.adj, p.signif.fdr)) |>
    left_join(y_position)

boxplot_ncells <- pd |>
    ggplot() +
    geom_boxplot(aes(y = ncells, x = APOE_carrier, fill = APOE_carrier), outlier.shape = NA) +
    geom_point(aes(y = ncells, x = APOE_carrier, fill = APOE_carrier)) +
    geom_text(data = n_cell_t_test, aes(label = fdr_anno, x = 1, y = y.position))+
    facet_wrap(~cell_type_anno, scales = "free_y") +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() 

ggsave(boxplot_ncells, filename = here(plot_dir, "boxplot_ncells_APOE_carrier.png"), height = 8 , width = 10)

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
# form <- ~ APOE + sample_id + exp_round
# form <- ~ APOE + Sex + Age + Anc_Afr + exp_round + seq_round + ncells

## TODO test  + pseudo_sum_umi +pseudo_expr_chrM, pseudo_expr_chrM_ratio

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