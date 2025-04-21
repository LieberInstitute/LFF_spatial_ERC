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
plot_data <- as.data.frame(colData(spe_pb)) |>
    select(sample_id, SpD, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, Visium_slide, round, nspots = ncells, pseudo_sum_umi, pseudo_expr_chrM_ratio) |>
    mutate(APOE = factor(APOE))
           

# Generate n sample boxplots for each SpD

sample_count <- plot_data |> group_by(APOE_carrier, SpD) |> count()

pb_SpD_bar_APOE_carrier <- plot_data |>
    count(APOE_carrier, SpD) |> 
    ggplot(aes(x = SpD, y=n, fill = APOE_carrier)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_carrier_colors)

ggsave(pb_SpD_bar_APOE_carrier, filename = here(plot_dir, 'pb_SpD_bar_APOE_carrier.png'))

pb_SpD_bar_APOE <- plot_data |>
    count(APOE, SpD) |> 
    ggplot(aes(x = SpD, y=n, fill = APOE)) +
    geom_col() +
    geom_text(aes(label = n), position = position_stack(vjust = .5)) +
    coord_flip() +
    scale_fill_manual(values = APOE_genotype_colors)

ggsave(pb_SpD_bar_APOE, filename = here(plot_dir, 'pb_SpD_bar_APOE.png'))

# Generate n cell boxplots for each SpD

pdf(here(plot_dir, "ERC_Visium_nspots_vs_APOE.pdf"), width = 8, height = 6)
for (d in levels(plot_data$SpD)) {
    
    SpD_data <- plot_data |> filter(SpD == d)
    y_max_ncells <- max(SpD_data$nspots)*1.1
    
    plot <- ggboxplot(
        SpD_data, 
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

form <- ~ APOE + sample_id + Sex + Age + Anc_Afr + Rin + Visium_slide + round + ncells

## TODO test  + pseudo_sum_umi +pseudo_expr_chrM, pseudo_expr_chrM_ratio

C <- canCorPairs(form, pd_sample)
# APOE and sample_id
# sample_id and Sex
# sample_id and Age
# sample_id and Anc_Afr
# sample_id and Rin
# sample_id and Visium_slide
# sample_id and round
# sample_id and ncells
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