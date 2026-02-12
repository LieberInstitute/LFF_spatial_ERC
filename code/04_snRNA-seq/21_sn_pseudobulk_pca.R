# Louise Huuki-Myers, March 2025
# Calculate and plot PCs for pseudobulked snRNA-seq data

library("SpatialExperiment")
library("scater")
library("tidyverse")
library("GGally")
library("here")
library("viridis")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("name", "n", "1", "character", "Data nickname",
      "path", "p", "2", "character", "path_to_data"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# opt$name <- "DE_input"
# opt$path <- "08_pseudoBulkDGE_sn/01_pseudobulk_data_sn/sce_pseudo_DGE-cell_type_anno.RDS"
# opt$cell_type = "Oligo.3"

plot_dir <- here("plots", "04_snRNA-seq", "21_sn_pseudobulk_pca", opt$name)
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# data_dir <- here("processed-data", "04_snRNA-seq", "21_sn_pseudobulk_pca")
# if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

## load colors
# load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### load pseuodbulked data ####
sce_pb <- readRDS(here("processed-data", opt$path))
cell_type_colors <- metadata(sce_pb)$cell_type_colors

## Drop sample Br1289
sce_pb <- sce_pb[,sce_pb$BrNum != 'Br1289']

#### Run PCA ####
set.seed(250320)

sce_pb <- scater::runPCA(sce_pb, 
                         ncomponents = 10,
                         name = "PCA")

pca <- reducedDim(sce_pb, "PCA")

## percent var
percent_var_df <- data.frame(PC = seq(ncol(pca)),
                             precentVar = attr(pca, "percentVar"))

pca_elbow <- ggplot(percent_var_df, aes(PC, precentVar)) +
    geom_line() +
    theme_bw()

ggsave(pca_elbow, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-PCA_elbow.png", opt$name)))

## add phenotype info to pca table

pd_select <- colData(sce_pb) |>
    as.data.frame() |>
    select(sample_id, APOE, Sex, Age, Ancestry, Anc_Afr, cell_type_anno, ncells, exp_round, seq_round) |>
    rownames_to_column("pseudobulk_sample")

pd_long_num <- pd_select |>
    select(pseudobulk_sample, Age,  Anc_Afr,  ncells) |>
    pivot_longer(!pseudobulk_sample, names_to = "var_num", values_to = "value_num")

pd_long_cat <- pd_select |>
    select(pseudobulk_sample, APOE, Sex, Ancestry, cell_type_anno, exp_round, seq_round) |>
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

ggsave(pca_vs_num_cor, filename = here(plot_dir , sprintf("sn_pseudobulk-%s-pca_vs_num_cor.png", opt$name)), height = 4, width = 4)

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

ggsave(pca_vs_cat_aov, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-pca_vs_cat_aov.png", opt$name)))


pca_vs_cat_boxplot <- pca_long |>
    ggplot(aes(x = value_cat, y = pca_value)) +
    geom_boxplot() +
    facet_grid(PCA~var_cat, scales = "free", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(pca_vs_cat_boxplot, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-pca_vs_cat_boxplot.png", opt$name)), height = 10, width = 10)

#### ggpairs ####
## TODO fix 
# Error in cor.test.default(x, y, method = method, use = use) : 
#     not enough finite observations

walk2(c("cell_type_anno", "APOE", "Sex", "Ancestry"), list(cell_type_colors$anno, APOE_genotype_colors, sex_colors, ancestry_colors),
      function(var, colors){
          
          gg_pca_plot <- ggpairs(pca_wide, columns = paste0("PC", 1:5), aes(colour = !!sym(var))) +
              scale_color_manual(values = colors) +
              scale_fill_manual(values = colors) +
              theme_bw()
          
          ggsave(gg_pca_plot, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-pca_ggpair_%s.png", opt$name, var)), height = 10, width = 10)
          
      })


#### select scatter plots ####

PC1_PC2_SpD <- pca_wide |>
    ggplot(aes(x = PC1, y = PC2, color = cell_type_anno, shape = APOE)) +
    geom_point() +
    scale_color_manual(values = cell_type_colors$anno) +
    theme_bw()

ggsave(PC1_PC2_SpD, filename = here(plot_dir , sprintf("sn_pseudobulk-%s-PC1_vs_PC2.png", opt$name)), width = 10)

PC1_PC6_SpD <- pca_wide |>
    ggplot(aes(x = PC1, y = PC6, color = cell_type_anno, shape = Sex)) +
    geom_point() +
    scale_color_manual(values = cell_type_colors$anno) +
    theme_bw()

ggsave(PC1_PC6_SpD, filename = here(plot_dir , sprintf("sn_pseudobulk-%s-PC1_vs_PC6_sex.png", opt$name)))


#### PCA by cell type ####

select_ct <- "Oligo.3"
sce_pb_select <- sce_pb[,sce_pb$cell_type_anno == select_ct]

sce_pb_select$cell_type_anno <- droplevels(sce_pb_select$cell_type_anno)

sce_pb_select <- scater::runPCA(sce_pb_select, 
                         ncomponents = 10,
                         name = "PCA")

pca <- reducedDim(sce_pb_select, "PCA")

dim(pca)

## percent var
percent_var_df <- data.frame(PC = seq(ncol(pca)),
                             precentVar = attr(pca, "percentVar"))

pca_elbow <- ggplot(percent_var_df, aes(PC, precentVar)) +
    geom_line() +
    theme_bw()

ggsave(pca_elbow, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-%s-PCA_elbow.png", opt$name, select_ct)))

## add phenotype info to pca table

pd_select <- colData(sce_pb_select) |>
    as.data.frame() |>
    select(sample_id, APOE, APOE_carrier, Sex, Age, Ancestry, Anc_Afr, cell_type_anno, ncells, exp_round, seq_round) |>
    rownames_to_column("pseudobulk_sample")

samples_of_interest <- c("Br5529","Br6423", "Br3974", "Br6263")

pca_wide <- cbind(pd_select, pca)  |>
    mutate(sample_oi = sample_id %in% samples_of_interest)

## explore cor

pd_long_num <- pd_select |>
    select(pseudobulk_sample, Age,  Anc_Afr,  ncells) |>
    pivot_longer(!pseudobulk_sample, names_to = "var_num", values_to = "value_num")

pd_long_cat <- pd_select |>
    select(pseudobulk_sample, APOE, Sex, Ancestry, cell_type_anno, exp_round, seq_round) |>
    pivot_longer(!pseudobulk_sample, names_to = "var_cat", values_to = "value_cat")

pca_long <- reshape2::melt(pca) |>
    rename(pseudobulk_sample = Var1, PCA = Var2, pca_value = value) |>
    left_join(pd_long_num, relationship = "many-to-many")|>
    left_join(pd_long_cat, relationship = "many-to-many")

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

ggsave(pca_vs_num_cor, filename = here(plot_dir , sprintf("sn_pseudobulk-%s-pca_vs_num_cor.png", opt$name)), height = 4, width = 4)


## select scatter 

library(ggrepel)

PC1_PC2_scatter <- pca_wide |>
    ggplot(aes(x = PC1, y = PC2, color = APOE, shape = Sex)) +
    geom_point() +
    geom_text_repel(aes(label = sample_id, color = sample_oi)) +
    scale_color_manual(values = c(APOE_genotype_colors, `TRUE` = "red") ) +
    theme_bw()

ggsave(PC1_PC2_scatter, filename = here(plot_dir , sprintf("sn_pseudobulk-%s-%s-PC1_vs_PC2.png", opt$name, select_ct)))

PC1_PC3_scatter <- pca_wide |>
    ggplot(aes(x = PC1, y = PC3, color = APOE_carrier, shape = Sex, size = ncells)) +
    geom_point() +
    geom_text_repel(aes(label = sample_id, color = sample_oi), size = 2) +
    scale_color_manual(values = c(APOE_carrier_colors, `TRUE` = "red") ) +
    theme_bw()

ggsave(PC1_PC3_scatter, filename = here(plot_dir , sprintf("sn_pseudobulk-%s-%s-PC1_vs_PC3.png", opt$name, select_ct)))



walk2(c("APOE", "Sex", "Ancestry", "sample_oi"), list(APOE_genotype_colors, sex_colors, ancestry_colors, c(`TRUE` = "red", `FALSE` = "grey20")),
      function(var, colors){
          
          gg_pca_plot <- ggpairs(pca_wide, columns = paste0("PC", 1:5), aes(colour = !!sym(var))) +
              scale_color_manual(values = colors) +
              scale_fill_manual(values = colors) +
              theme_bw()
          
          ggsave(gg_pca_plot, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-%s-pca_ggpair_%s.png", opt$name, select_ct, var)), height = 10, width = 10)
          
      })

# walk(c("Age", "ncells", "Anc_Afr"),
#       function(var){
#           
#           gg_pca_plot <- ggpairs(pca_wide, columns = paste0("PC", 1:5), aes(color = Age)) +
#               # scale_color_manual(values = colors) +
#               # scale_fill_manual(values = colors) +
#               theme_bw()
#           
#           ggsave(gg_pca_plot, filename = here(plot_dir, sprintf("sn_pseudobulk-%s-%s-pca_ggpair_%s.png", opt$name, select_ct, var)), height = 10, width = 10)
#           
#       })


# slurmjobs::job_single('21_sn_pseudobulk_pca', create_shell = TRUE, memory = '10G', command = "Rscript 21_sn_pseudobulk_pca.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
