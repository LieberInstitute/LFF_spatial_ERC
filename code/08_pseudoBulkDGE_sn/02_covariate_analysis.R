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
library("variancePartition")
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

table(sce_pb$cell_type_anno)

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

colnames(colData(sce_pb))

# Extract the relevant columns from the data
pd <- as.data.frame(colData(sce_pb))

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

## Drop cell types with low samples (<3)
(low_n_ct <- cell_type_count |> filter(n < 3) |> pull(cell_type_anno) |> unique())
# [1] Excit.L5_6_NP Inhib.Pax6  

sce_pb <- sce_pb[,! sce_pb$cell_type_anno %in% low_n_ct]
sce_pb$cell_type_anno <- droplevels(sce_pb$cell_type_anno)

table(sce_pb$cell_type_anno)
pd <- as.data.frame(colData(sce_pb))

#### t-test variables ####

test_variables <- c("ncells", "pseudo_sum_umi", "pseudo_expr_chrM_ratio")
names(test_variables) <- test_variables

summary(pd[,test_variables])

var_t_test <- map(test_variables, function(test_var){
    
    y_position <- pd |>
        group_by(cell_type_anno) |>
        summarise(y.position = max(!!sym(test_var)) + .05*max(!!sym(test_var)))
    
    var_t_test <- pd |>
        do(compare_means(!!sym(test_var) ~ APOE_carrier, data = ., method = "t.test", p.adjust.method = "fdr", group.by = "cell_type_anno")) |>
        ungroup() |>
        mutate(p.signif.fdr = case_when(p.adj < 0.005 ~ "***",
                                        p.adj < 0.01 ~"**",
                                        p.adj < 0.05 ~"*",
                                        TRUE~""),
               fdr_anno = sprintf("FDR=%.3f%s", p.adj, p.signif.fdr)) |>
        left_join(y_position)
    
    boxplot_test_var <- pd |>
        ggplot() +
        geom_boxplot(aes(y = !!sym(test_var), x = APOE_carrier, fill = APOE_carrier), outlier.shape = NA) +
        geom_point(aes(y = !!sym(test_var), x = APOE_carrier, fill = APOE_carrier)) +
        geom_text(data = var_t_test, aes(label = fdr_anno, x = 1, y = y.position))+
        facet_wrap(~cell_type_anno, scales = "free_y") +
        scale_fill_manual(values = APOE_carrier_colors) +
        theme_bw() 
    
    ggsave(boxplot_test_var, filename = here(plot_dir, sprintf("boxplot_%s_APOE_carrier.png", test_var)), width = 10)
    
    return(var_t_test)
    
})

map(var_t_test, ~.x|> filter(p.adj < 0.05))


#### variance partition on Sample level data ####

# load sample level data
dim(sce_pb)
table(sce_pb$APOE)
# E2/E2 E2/E3 E3/E4 E4/E4 
# 86   142   173   123
table(sce_pb$APOE_carrier)
# E2+ E4+ 
# 228 296

table(is.na(pd$Rin))

# Assess correlation between all pairs of variables
colnames(pd)

form <- ~ cell_type_anno + APOE + sample_id + Sex + Age + Anc_Afr + Rin + exp_round + seq_round + ncells + pseudo_sum_umi + pseudo_expr_chrM_ratio
# form <- ~ APOE + sample_id + exp_round
# form <- ~ APOE + Sex + Age + Anc_Afr + exp_round + seq_round + ncells

C <- canCorPairs(form, pd)
# ~ APOE + sample_id + Sex + Age + Anc_Afr + exp_round + seq_round + ncells
# APOE and sample_id
# sample_id and Sex
# sample_id and Age
# sample_id and Anc_Afr
# sample_id and Rin
# sample_id and exp_round
# sample_id and seq_round
# sample_id and ncells

pdf(here(plot_dir, "sn_variable_cor_matrix.pdf"))
plotCorrMatrix(C)
dev.off()

## test variance partition 
## TODO add Rin

message(Sys.time(), " - fitExtractVarPartModel")

# form <- ~ (1 | cell_type_anno) + (1 | APOE_carrier) + (1 | APOE) + (1 | Sex) + Anc_Afr + Age + (1 | exp_round) + (1 | seq_round) + pseudo_sum_umi + pseudo_expr_chrM_ratio
form <- ~ (1 | cell_type_anno) + (1 | APOE_carrier) + (1 | APOE) + (1 | Sex) + Anc_Afr + Age + (1 | exp_round) + (1 | seq_round)
varPart <- fitExtractVarPartModel(logcounts(sce_pb), form, pd)
vp <- sortCols(varPart)

## violin plot 
sn_vp_violin <- plotVarPart(vp)
ggsave(sn_vp_violin, filename = here(plot_dir, "sn_vp_violin.png"))

# slurmjobs::job_single('02_covariate_analysis', create_shell = TRUE, memory = '50G', command = "Rscript 02_covariate_analysis.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()