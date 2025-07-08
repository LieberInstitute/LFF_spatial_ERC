## Louise Huuki-Myers, July 2025
## Multicellular factor analysis
## adapted from https://github.com/LieberInstitute/MFA_LC-ERC_pilot/blob/a99c526d6dfa1f97b274b4d3f65de7da75f0c299/code/01_main_analysis/01_mofa.R

#### set up ####
library("SpatialExperiment")
library("MOFAcellulaR")
library("MOFA2")
library("tidyverse")
library("GGally")
library("here")
library("sessioninfo")

# out_path = here('processed-data', '01_main_analysis', 'model.hdf5')
# plot_dir = here('plots', '01_main_analysis', 'MOFA')

data_dir <- here("processed-data", "14_MOFA", "01_MOFA")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "14_MOFA", "01_MOFA")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## colors
load(here("processed-data", "project_colors.Rdata"))
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)


set.seed(1)

#### Load & Combine Data ####

# erc_cluster_var = 'BayesSpace_PCA_Harmony_k09'
# lc_cluster_var = 'cluster_HARMONYlmbna_hvg20k10'

#   Load SPE objects
visium_spe <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))
sn_sce <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS"))

#   Subset to shared genes
shared_genes <- intersect(rownames(sn_sce), rownames(visium_spe))
visium_spe <- visium_spe[shared_genes, ]
sn_sce <- sn_sce[shared_genes, ]

nrow(visium_spe) == nrow(sn_sce)

# Add data type
visium_spe$data_type <- "Visium"
sn_sce$data_type <- "snRNA-seq"

## create new sample ID cols
visium_spe$sample_cluster <- paste0(visium_spe$BrNum, "_", visium_spe$registration_variable)
sn_sce$sample_cluster <- paste0(sn_sce$BrNum, "_", sn_sce$registration_variable)

## update colnames
colnames(visium_spe) <- visium_spe$sample_cluster 
colnames(sn_sce) <- sn_sce$sample_cluster 

visium_spe$sample_id <- paste0(visium_spe$sample_id, "_Visium")
sn_sce$sample_id <- paste0(sn_sce$sample_id, "_snRNA-seq")

common_col_data <- intersect(colnames(colData(sn_sce)),colnames(colData(visium_spe)))

colData(sn_sce) <- colData(sn_sce)[,common_col_data]
colData(visium_spe) <- colData(visium_spe)[,common_col_data]

## drop reduced dims and metadata
reducedDims(sn_sce) = list()
reducedDims(visium_spe) = list()

metadata(sn_sce) = list()
metadata(visium_spe) = list()

## check colData
identical(colnames(colData(sn_sce)),colnames(colData(visium_spe)))
## match factors
colData(visium_spe)$APOE <- factor(colData(visium_spe)$APOE, levels = levels(colData(sn_sce)$APOE))
colData(visium_spe)$APOE_carrier <- factor(colData(visium_spe)$APOE_carrier, levels = levels(colData(sn_sce)$APOE_carrier))
colData(visium_spe)$Ancestry <- factor(colData(visium_spe)$Ancestry, levels = levels(colData(sn_sce)$Ancestry))
colData(visium_spe)$Sex <- factor(colData(visium_spe)$Sex, levels = levels(colData(sn_sce)$Sex))

v_class <- map_chr(colData(visium_spe), class)
s_class <- map_chr(colData(sn_sce), class)

all(v_class == s_class)

## sync rowData
rowData(sn_sce) <- rowData(visium_spe)
rowRanges(sn_sce) <- rowRanges(visium_spe)

assayNames(sn_sce)
assayNames(visium_spe)

## convert spe to sce
visium_spe <- as(visium_spe, "SingleCellExperiment")
class(sn_sce)
class(visium_spe)

# spe = cbind(visium_spe, sn_sce)
# Error in value[[3L]](cond) : 
#     failed to combine 'int_colData' in 'cbind(<SingleCellExperiment>)':
#     the DFrame objects to combine must have the same column names

spe <- SingleCellExperiment(colData = rbind(colData(visium_spe), colData(sn_sce)),
                            rowData = rowData(visium_spe),
                            assays = list(logcounts = cbind(logcounts(visium_spe),
                                                            logcounts(sn_sce)),
                                          counts = cbind(counts(visium_spe),
                                                            counts(sn_sce)))
                            )

spe$cluster <- spe$registration_variable

## get donor data
pd <- colData(spe) |>
    as_tibble() |>
    select(BrNum, Age, Sex, APOE, Ancestry, APOE_carrier) |>
    unique()

#### Preprocess expression and create a MOFA object ####
set.seed(708)

mofa = create_init_exp(
    counts = assays(spe)$counts, coldata = as.data.frame(colData(spe))
) |>
    filt_profiles(
        #   Retain all spatial domains initially
        cts = unique(spe$cluster),
        ncells = 0,
        counts_col = "ncells",
        ct_col = "cluster"
    ) |>
    #   Drop spatial domains seen in very few samples
    filt_views_bysamples(nsamples = 2) |>
    #   Require a gene to have 5 counts in at least 25% of samples
    filt_gex_byexpr(min.count = 5, min.prop = 0.25) |>
    #   Drop spatial domains having below 15 genes at this point
    filt_views_bygenes(ngenes = 15) |>
    #   Drop samples below 90% coverage of genes
    filt_samples_bycov(prop_coverage = 0.9) |>
    #   Normalize expression
    tmm_trns(scale_factor = 1000000) |>
    #   Filter to HVGs only
    filt_gex_byhvg(prior_hvg = NULL, var.threshold = 0) |>
    #   Again, drop spatial domains having below 15 genes at this point
    filt_views_bygenes(ngenes = 15) |>
    pb_dat2MOFA(sample_column = 'BrNum') |>
    create_mofa()

#   Metadata seems to need to be manually added  
samples_metadata(mofa) = left_join(
    samples_metadata(mofa), pd, by = c("sample" = "BrNum")
)

####   Fit the MOFA model ####

#   Use mostly default options except where overriden in
#   https://saezlab.github.io/MOFAcellulaR/articles/get-started.html#fitting-a-mofa-model.
#   Use 5 factors in response to warnings when using any more
data_opts = get_default_data_options(mofa)
train_opts = get_default_training_options(mofa)
model_opts = get_default_model_options(mofa)
model_opts$spikeslab_weights = FALSE 
model_opts$num_factors = 5

mofa <- prepare_mofa(
    object = mofa,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

model <- run_mofa(mofa, data_dir, use_basilisk = TRUE)


####  Exploratory plots ####
#   Weights to each factor grouped by APOE genotype
p = plot_factor(
    model, factors = "all", color_by = 'APOE_geno', dodge = TRUE,
    add_violin = TRUE
)
pdf(file.path(plot_dir, 'weights_by_APOE.pdf'), width = 10, height = 5)
print(p)
dev.off()

#   Get factor weights for each donor
factor_df = get_tidy_factors(
    model = model,
    metadata = samples_metadata(model),
    factor = "all",
    sample_id_column = "sample"
)

#   For each factor, compute p-value of the coefficients in linear regression
#   of weight value against APOE numeric genotype
p_val_df_list = list()
for (this_factor in paste0("Factor", 1:5)) {
    this_factor_df = factor_df |>
        filter(Factor == this_factor)
    
    this_lm = lm(value ~ APOE, data = this_factor_df)
    p_val_df_list[[this_factor]] = tibble(
        Factor = this_factor,
        p_value = signif(
            summary(this_lm)$coefficients["APOE", 4], 2
        )
    )
}
p_val_df = do.call(rbind, p_val_df_list)

#   Manually determine where to place p-value labels vertically
p_val_df$y = c(Inf, -Inf, Inf, -Inf, -Inf)
p_val_df$vjust = c(1, 0, 1, 0, 0)

#   Weights to each factor grouped by APOE genotype (figure-ready version)
p = ggplot(factor_df, mapping = aes(x = APOE_geno, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_smooth(mapping = aes(group = 1), method = "lm", se = FALSE) +
    geom_jitter(size = 3, width = 0.1) +
    geom_text(
        data = p_val_df,
        mapping = aes(
            x = Inf, y = y, label = sprintf("\np = %s  \n", p_value),
            vjust = vjust
        ),
        hjust = 1, size = 10
    ) +
    facet_wrap(~Factor, nrow = 1, scales = "free") +
    theme_bw(base_size = 30) +
    theme(
        axis.text.x = element_text(
            angle = 90, vjust = 0.5, hjust = 1
        )
    ) +
    labs(x = "APOE Genotype", y = "Weight")
pdf(file.path(plot_dir, 'weights_by_APOE_boxplot.pdf'), width = 28, height = 7)
print(p)
dev.off()

#   Convert to wide format
factor_df = factor_df |>
    pivot_wider(names_from = "Factor", values_from = "value") |>
    relocate(matches('^Factor'))

#   Relationship between factors colored by presence of E4 allele
p = ggpairs(factor_df, columns = 1:5, aes(color = factor(APOE_E4_pos)))
pdf(file.path(plot_dir, 'ggpairs_APOE_E4_pos.pdf'))
print(p)
dev.off()

#   Relationship between factors colored by APOE genotype
p = ggpairs(factor_df, columns = 1:5, aes(color = factor(APOE_geno)))
pdf(file.path(plot_dir, 'ggpairs_APOE_geno.pdf'))
print(p)
dev.off()

#   Test association of APOE genotype with each factor
assoc_list = list()
assoc_list[['APOE']] = get_associations(
    model = model,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    test_variable = "APOE",
    test_type = "continuous",
    group = FALSE
)
print(assoc_list)

#   Show unadjusted p-value for consistency with other plots
assoc_list[['APOE']]$adj_pvalue = assoc_list[['APOE']]$p.value

#   Plot a heatmap of summary results, labeling with covariates of interest
p = plot_MOFA_hmap(
    model = model,
    group = FALSE,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    sample_anns = c("APOE_geno", "Ancestry", "Age", "Sex"),
    assoc_list = assoc_list,
    col_rows = list(
        'APOE_geno' = APOE_geno_colors,
        'Ancestry' = ancestry_colors,
        'Sex' = sex_colors
    )
)
pdf(file.path(plot_dir, 'heatmap.pdf'))
print(p)
dev.off()

session_info()

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/