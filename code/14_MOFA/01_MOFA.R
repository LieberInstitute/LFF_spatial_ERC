library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(GGally)

out_path = here('processed-data', '01_main_analysis', 'model.hdf5')
plot_dir = here('plots', '01_main_analysis', 'MOFA')

erc_spe_path = here(
    'processed-data', '00_project_prep', '01_spe_model_pseudobulk',
    'pseudobulk', 'spe_ERC_pseudobulk-spatial_cluster.rds'
)
erc_cluster_var = 'BayesSpace_PCA_Harmony_k09'
lc_spe_path = here(
    'processed-data', '00_project_prep', '01_spe_model_pseudobulk',
    'pseudobulk', 'spe_LC_pseudobulk-spatial_cluster.rds'
)
lc_cluster_var = 'cluster_HARMONYlmbna_hvg20k10'

#   Custom colors for donor-level metadata
APOE_geno_colors <- c(
    `E2/E2`="#114B5F", `E2/E3` = "#61C9A8", `E3/E4`="#ED9B40", `E4/E4`="#BA3B46"
)
ancestry_colors <- c(EA = "#1B3174", AA = "#698F3F")
sex_colors <- c(M = "#5C80BC", F = "#D58BCC")

dir.create(plot_dir, showWarnings = FALSE)
set.seed(1)

################################################################################
#   Combine LC and ERC SPE objects
################################################################################

#   Load SPE objects
lc_spe = readRDS(lc_spe_path)
erc_spe = readRDS(erc_spe_path)

#   Subset to shared genes
shared_genes = intersect(rownames(erc_spe), rownames(lc_spe))
lc_spe = lc_spe[shared_genes, ]
erc_spe = erc_spe[shared_genes, ]

#   Gather some donor-level info for later
pd = erc_spe[, match(unique(erc_spe$BrNum), erc_spe$BrNum)] |>
    colData() |>
    as_tibble() |>
    select(BrNum, Age, Sex, APOE, Ancestry) |>
    rename(APOE_geno = APOE) |>
    mutate(
        APOE_geno = factor(
            sub('\\.', '/', APOE_geno),
            levels = c("E2/E2", "E2/E3", "E3/E4", "E4/E4")
        ),
        APOE = as.integer(APOE_geno),
        APOE_E4_pos = grepl('E4', APOE_geno)
    )

#   Simplify ERC SPE object in preparation for combining with LC
erc_spe$cluster = sub('Sp09D0', 'ERC_D', erc_spe[[erc_cluster_var]])
erc_spe$brain_region = 'ERC'
colData(erc_spe) = colData(erc_spe) |>
    as_tibble() |>
    select(sample_id, BrNum, brain_region, cluster, APOE, ncells) |>
    mutate(sample_id_better = paste(BrNum, cluster, sep = '_')) |>
    DataFrame()
colnames(erc_spe) = erc_spe$sample_id_better
reducedDims(erc_spe) = list()

#   Simplify LC SPE object in preparation for combining with ERC
lc_spe$cluster = paste0('LC_D', lc_spe[[lc_cluster_var]])
lc_spe$brain_region = 'LC'
colData(lc_spe) = colData(lc_spe) |>
    as_tibble() |>
    select(sample_id, BrNum, brain_region, cluster, APOE, ncells) |>
    mutate(sample_id_better = paste(BrNum, cluster, sep = '_')) |>
    DataFrame()
colnames(lc_spe) = lc_spe$sample_id_better
reducedDims(lc_spe) = list()

#   Merge regions, with the 'cluster' variable now including info about region
#   and spatial domain
spe = cbind(lc_spe, erc_spe)

################################################################################
#   Preprocess expression and create a MOFA object
################################################################################

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

################################################################################
#   Fit the MOFA model
################################################################################

#   Use mostly default options except where overriden in
#   https://saezlab.github.io/MOFAcellulaR/articles/get-started.html#fitting-a-mofa-model.
#   Use 5 factors in response to warnings when using any more
data_opts = get_default_data_options(mofa)
train_opts = get_default_training_options(mofa)
model_opts = get_default_model_options(mofa)
model_opts$spikeslab_weights = FALSE 
model_opts$num_factors = 5

mofa = prepare_mofa(
    object = mofa,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

model = run_mofa(mofa, out_path, use_basilisk = TRUE)

################################################################################
#    Exploratory plots
################################################################################

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