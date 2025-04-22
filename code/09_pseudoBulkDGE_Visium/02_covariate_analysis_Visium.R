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
library("variancePartition")

#### Set up dirs ####
plot_dir <- here("plots", "09_pseudoBulkDGE_Visium", "02_covariate_analysis_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load the data ####
spe_pb <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))

## load colors
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

# Extract the relevant columns from the data
pd <- as.data.frame(colData(spe_pb)) |>
    select(sample_id, SpD, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, Visium_slide, round, nspots = ncells, pseudo_sum_umi, pseudo_expr_chrM_ratio) |>
    mutate(APOE = factor(APOE))
           

#### t-test variables ####

## nspots
y_position <- pd |>
    group_by(SpD) |>
    summarise(y.position = max(nspots) + .05*max(nspots))

nspots_t_test <- pd |>
    do(compare_means(nspots ~ APOE_carrier, data = ., method = "t.test", p.adjust.method = "fdr", group.by = "SpD")) |>
    ungroup() |>
    mutate(p.signif.fdr = case_when(p.adj < 0.005 ~ "***",
                                    p.adj < 0.01 ~"**",
                                    p.adj < 0.05 ~"*",
                                    TRUE~""),
           fdr_anno = sprintf("FDR=%.3f%s", p.adj, p.signif.fdr)) |>
    left_join(y_position)

boxplot_nspots <- pd |>
    ggplot() +
    geom_boxplot(aes(y = nspots, x = APOE_carrier, fill = APOE_carrier), outlier.shape = NA) +
    geom_point(aes(y = nspots, x = APOE_carrier, fill = APOE_carrier)) +
    geom_text(data = nspots_t_test, aes(label = fdr_anno, x = 1, y = 1500))+
    facet_wrap(~SpD) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() 

ggsave(boxplot_nspots, filename = here(plot_dir, "boxplot_nspots_APOE_carrier.png"), height = 8 , width = 10)

## other variables
t_test_variables <- c("pseudo_sum_umi", "pseudo_expr_chrM_ratio")
names(t_test_variables) <- t_test_variables
var_t_test <- map(t_test_variables,
                  ~pd |>
                      do(compare_means(!!sym(.x) ~ APOE_carrier, data = ., method = "t.test", p.adjust.method = "fdr", group.by = "SpD")) |>
                      mutate(p.signif.fdr = case_when(p.adj < 0.005 ~ "***",
                                                      p.adj < 0.01 ~"**",
                                                      p.adj < 0.05 ~"*",
                                                      TRUE~""),
                             fdr_anno = sprintf("FDR=%.3f%s", p.adj, p.signif.fdr)) 
                  )





#### variance parition on Sample level data ####

spe_pb_sample <- readRDS(here("processed-data", "05_spe_correct_cluster", "25_SpD_pseudobulk", "spe_pseudobulk_only-sample_id.rds"))
dim(spe_pb_sample)
table(spe_pb_sample$APOE)
# E2/E2 E2/E3 E3/E4 E4/E4 
# 6     8    10     7 
table(spe_pb_sample$APOE_carrier)
# E2+ E4+ 
# 14  17 

pd_sample <- as.data.frame(colData(spe_pb_sample))

# Assess correlation between all pairs of variables
colnames(pd_sample)

form <- ~ APOE + sample_id + Sex + Age + Anc_Afr + Rin + Visium_slide + round + nspots

## TODO test  + pseudo_sum_umi +pseudo_expr_chrM, pseudo_expr_chrM_ratio

C <- canCorPairs(form, pd_sample)
# APOE and sample_id
# sample_id and Sex
# sample_id and Age
# sample_id and Anc_Afr
# sample_id and Rin
# sample_id and Visium_slide
# sample_id and round
# sample_id and nspots
# Visium_slide and round

pdf(here(plot_dir, "Visium_variable_cor_matrix.pdf"))
plotCorrMatrix(C)
dev.off()

## test variance partition 
form <- ~ (1 | APOE_carrier) + (1 | APOE) + (1 | Sex) + Anc_Afr + Age + Rin + (1 | Visium_slide) + (1 | round)
varPart <- fitExtractVarPartModel(logcounts(spe_pb_sample), form, pd_sample)
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