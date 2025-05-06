# Louise Huuki-Myers, March 2025
# Calculate and plot PCs for pseudobulked Visium data

library("SpatialExperiment")
library("scater")
library("tidyverse")
library("GGally")
library("here")
library("viridis")
library("sessioninfo")

plot_dir <- here("plots", "05_spe_correct_cluster", "28_SpD_pseudobulk_pca")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "05_spe_correct_cluster", "28_SpD_pseudobulk_pca")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load colors
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### load pseuodbulked data ####
spe_pb <- readRDS(here("processed-data","05_spe_correct_cluster", "20_model_pseudobulk_anno","spe_pseudobulk-SpD.rds"))


#### Run PCA ####
set.seed(250320)

spe_pb <- scater::runPCA(spe_pb, 
                         ncomponents = 10,
                         name = "PCA")

pca <- reducedDim(spe_pb, "PCA")

## percent var
percent_var_df <- data.frame(PC = seq(ncol(pca)),
                             precentVar = attr(pca, "percentVar"))

pca_elbow <- ggplot(percent_var_df, aes(PC, precentVar)) +
    geom_line() +
    theme_bw()

ggsave(pca_elbow, filename = here(plot_dir, "PCA_elbow.png"))

## add phenotype info to pca table

pd_select <- colData(spe_pb) |>
    as.data.frame() |>
    select(sample_id, APOE, Sex, Age, Ancestry, Anc_Afr, SpD, ncells, round, Visium_slide) |>
    rownames_to_column("pseudobulk_sample")

pd_long_num <- pd_select |>
    select(pseudobulk_sample, Age,  Anc_Afr,  ncells) |>
    pivot_longer(!pseudobulk_sample, names_to = "var_num", values_to = "value_num")

pd_long_cat <- pd_select |>
    select(pseudobulk_sample, APOE, Sex, Ancestry, SpD, round, Visium_slide) |>
    pivot_longer(!pseudobulk_sample, names_to = "var_cat", values_to = "value_cat")

pca_long <- reshape2::melt(pca) |>
    rename(pseudobulk_sample = Var1, PCA = Var2, pca_value = value) |>
    left_join(pd_long_num, relationship = "many-to-many")|>
    left_join(pd_long_cat, relationship = "many-to-many")

pca_wide <- cbind(pd_select, pca)

#### correlation and Anova ####

pca_cor <- pca_long |>
    group_by(PCA, var_num) |>
    summarize(cor = cor(pca_value, value_num))

pca_cor |> arrange(-abs(cor))

## plot correlations
pca_vs_num_cor <- pca_cor |>
    ggplot(aes(x = var_num, y = PCA, fill = cor)) +
    geom_tile() +
    geom_text(aes(label = round(cor, 2)), size = 2) +
    theme_bw() +
    scale_fill_gradient2() 

ggsave(pca_vs_num_cor, filename = here(plot_dir , "Visium_peudobulk_pca_vs_num_cor.png"), height = 4, width = 4)

## anova on cat variables
df_aov <- pca_long |>
    group_by(PCA, var_cat) |>
    do(aov = broom::tidy(aov(pca_value ~ value_cat, data = .)))|>
    unnest(aov) |>
    filter(term == "value_cat") |>
    mutate(FDR = p.adjust(p.value, method = "fdr"))

df_aov |> count(FDR < 0.001)

## plot anova p-values
pca_vs_cat_aov <- df_aov |>
    ggplot(aes(x = var_cat, y = PCA, fill = -log10(p.value))) +
    geom_tile() +
    # geom_text(aes(label = round(cor, 2)), size = 2) +
    theme_bw() +
    scale_fill_viridis()

ggsave(pca_vs_cat_aov, filename = here(plot_dir , "Visium_pseudobulk_pca_vs_cat_aov.png"))


pca_vs_cat_boxplot <- pca_long |>
    ggplot(aes(x = value_cat, y = pca_value)) +
    geom_boxplot() +
    facet_grid(PCA~var_cat, scales = "free", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(pca_vs_cat_boxplot, filename = here(plot_dir , "Visium_pseudobulk_pca_vs_cat_boxplot.png"), height = 10)

#### ggpairs ####
## TODO fix 
# Error in cor.test.default(x, y, method = method, use = use) : 
#     not enough finite observations

walk2(c("SpD", "APOE", "Sex", "Ancestry"), list(SpD_colors, APOE_genotype_colors, sex_colors, ancestry_colors),
      function(var, colors){
          
          gg_pca_plot <- ggpairs(pca_wide, columns = paste0("PC", 1:5), aes(colour = !!sym(var))) +
              scale_color_manual(values = colors) +
              scale_fill_manual(values = colors) +
              theme_bw()
          
          ggsave(gg_pca_plot, filename = here(plot_dir, sprintf("Visium_pseudobulk_pca_ggpair_%s.png", var)))
          
      })


#### select scatter plots ####

PC1_PC2_SpD <- pca_wide |>
    ggplot(aes(x = PC1, y = PC2, color = SpD, shape = APOE)) +
    geom_point() +
    scale_color_manual(values = SpD_colors) +
    theme_bw()

ggsave(PC1_PC2_SpD, filename = here(plot_dir , "Visium_pseudobulk_PC1_vs_PC2.png"))

PC1_PC6_SpD <- pca_wide |>
    ggplot(aes(x = PC1, y = PC6, color = SpD, shape = Sex)) +
    geom_point() +
    scale_color_manual(values = SpD_colors) +
    theme_bw()

ggsave(PC1_PC6_SpD, filename = here(plot_dir , "Visium_pseudobulk_PC1_vs_PC6_sex.png"))


# slurmjobs::job_single('28_SpD_pseudobulk_pca', create_shell = TRUE, memory = '10G', command = "Rscript 28_SpD_pseudobulk_pca.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
