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
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test 
# opt$datatype = "sn_broad"

data_dir <- here("processed-data", "14_MOFA", "01_MOFA", opt$datatype)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "14_MOFA", "01_MOFA", opt$datatype)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## colors
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

#### Load & Combine Data ####

#   Load SPE objects
visium_spe <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))

if(opt$datatype == "sn_fine"){
    message("sn FINE")
    sn_sce <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_anno.RDS"))
} else if(opt$datatype == "sn_broad"){
    message("sn BROAD")
    sn_sce <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS"))
}

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

table(spe$cluster)

rm(visium_spe)
rm(sn_sce)

## load tau pathology data
path_tb <- read.csv(here("processed-data", "00_project_prep","05_pathology","sample_taupathy.csv")) |>
    column_to_rownames("BrNum")

spe$taupathy <- ifelse(path_tb[spe$BrNum,]$taupathy, "t+", "t-")
spe$CERAD <- path_tb[spe$BrNum,]$CERAD
spe$Braak <- path_tb[spe$BrNum,]$Braak

## Drop Br1289
spe <- spe[,spe$BrNum != "Br1289"]

table(spe$BrNum, spe$data_type)

## get donor data
pd <- colData(spe) |>
    as_tibble() |>
    select(BrNum, Age, Sex, Ancestry, Anc_Afr, APOE, APOE_carrier, taupathy, CERAD, Braak) |>
    unique()

#### Preprocess expression and create a MOFA object ####
message(Sys.time() , " - Create MOFA object")
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
    filt_views_bysamples(nsamples = 2) |> #   Drop spatial domains seen in very few samples
    filt_gex_byexpr(min.count = 5, min.prop = 0.25) |> #   Require a gene to have 5 counts in at least 25% of samples
    filt_views_bygenes(ngenes = 15) |> #   Drop spatial domains having below 15 genes at this point
    filt_samples_bycov(prop_coverage = 0.9) |> #   Drop samples below 90% coverage of genes
    tmm_trns(scale_factor = 1000000) |> #   Normalize expression
    filt_gex_byhvg(prior_hvg = NULL, var.threshold = 0) |> #   Filter to HVGs only
    filt_views_bygenes(ngenes = 15) |>  #   Again, drop spatial domains having below 15 genes at this point
    pb_dat2MOFA(sample_column = 'BrNum') |>
    create_mofa()

#   Metadata seems to need to be manually added  
samples_metadata(mofa) = left_join(
    samples_metadata(mofa), pd, by = c("sample" = "BrNum")
)

####   Fit the MOFA model ####
message(Sys.time() , " - Fit the MOFA model")
#   Use mostly default options except where overriden in
#   https://saezlab.github.io/MOFAcellulaR/articles/get-started.html#fitting-a-mofa-model.

data_opts = get_default_data_options(mofa)
train_opts = get_default_training_options(mofa)
model_opts = get_default_model_options(mofa)
model_opts$spikeslab_weights = FALSE 
model_opts$num_factors = 7 #   Use 7 factors in response to warnings when using any more

mofa <- prepare_mofa(
    object = mofa,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

out_path = here(data_dir, 'model.hdf5')

model <- run_mofa(mofa, out_path, use_basilisk = TRUE)

## how to read in mofa?
# library(rhdf5)
# h5read(name = out_path)

####  Exploratory plots ####

#   Get factor weights for each donor - export to csv
factor_df = get_tidy_factors(
    model = model,
    metadata = samples_metadata(model),
    factor = "all",
    sample_id_column = "sample"
)

## check min and max samples
factor_df |>
    group_by(Factor) |>
    slice_max(value)

factor_df |>
    group_by(Factor) |>
    slice_min(value)

write_csv(factor_df, file = here(data_dir, sprintf("MOFA_factor_df-%s.csv", opt$datatype)))

#   Test association of APOE genotype with each factor

test_vars <- c('APOE_carrier', 'Ancestry', "taupathy", "Sex", "Age","Braak","CERAD")
names(test_vars) <- test_vars
    
assoc_list <- map(test_vars, ~get_associations( model = model,
                                                metadata = samples_metadata(model),
                                                sample_id_column = "sample",
                                                test_variable = .x,
                                                test_type = "categorical",
                                                group = FALSE))

assoc_tb <- do.call("bind_rows", assoc_list) |>
    mutate(signif = case_when(adj_pvalue < 0.001 ~ "***",
                              adj_pvalue < 0.01 ~ "**",
                              adj_pvalue < 0.05 ~ "*",
                              TRUE ~ "")
    ) 

assoc_tb |> group_by(term) |> arrange(adj_pvalue) |> slice(1)


#   Show unadjusted p-value for consistency with other plots ??
# assoc_list[['APOE']]$adj_pvalue = assoc_list[['APOE']]$p.value


factor_boxplot <- function(var, fill_colors = NULL, text45 = TRUE, assoc_tb){
    
    weights_boxplot <- ggplot(factor_df) +
        geom_boxplot(aes(x = !!sym(var), y = value, fill = !!sym(var)), outlier.shape = NA) +
        geom_jitter(aes(x = !!sym(var), y = value, fill = !!sym(var)), width = 0.1) +
        facet_wrap(~Factor, nrow = 1, scales = "free") +
        labs(y = "Weight") +
        theme_bw() +
        theme(legend.position = "None") 
    
        weights_boxplot <- weights_boxplot + ggplot2::geom_label(
            data = assoc_tb |> filter(term == var), 
            ggplot2::aes(x = -Inf, y = -Inf, label = sprintf("pval=%.2e%s", adj_pvalue, signif)),
            alpha = 0.5,
            vjust = "inward", 
            hjust = "inward", 
            size = 2.5
        )

    if(!is.null(fill_colors)) weights_boxplot <- weights_boxplot + scale_fill_manual(values = fill_colors)
    if(text45) weights_boxplot <- weights_boxplot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(weights_boxplot, filename = here(plot_dir, sprintf("factor_weights_boxplot_%s.png", var)), height = 5, width = 10)
    
}

factor_boxplot(var = "APOE_carrier", fill_colors = APOE_carrier_colors, assoc_tb = assoc_tb)
factor_boxplot(var = "Ancestry", fill_colors = ancestry_colors, assoc_tb = assoc_tb)
factor_boxplot(var = "Sex", fill_colors = sex_colors, assoc_tb = assoc_tb)
factor_boxplot(var = "taupathy", assoc_tb = assoc_tb)
factor_boxplot(var = "Braak", assoc_tb = assoc_tb)


#   Convert to wide format
factor_df = factor_df |>
    pivot_wider(names_from = "Factor", values_from = "value") |>
    relocate(matches('^Factor'))

#   Relationship between factors and APOE carrier
p = ggpairs(factor_df, columns = 1:5, aes(color = APOE_carrier))

pdf(file.path(plot_dir, 'ggpairs_APOE_carrier.pdf'))
print(p)
dev.off()

#   Relationship between factors colored by APOE genotype
p = ggpairs(factor_df, columns = 1:5, aes(color = factor(APOE)))
pdf(file.path(plot_dir, 'ggpairs_APOE_geno.pdf'))
print(p)
dev.off()


#   Plot a heatmap of summary results, labeling with covariates of interest
pdf(here(plot_dir, sprintf('MOFA_heatmap_%s.pdf', opt$datatype)), height = 3 + (length(model@data)/4))
plot_MOFA_hmap(
    model = model,
    group = FALSE,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    sample_anns = c("APOE_carrier", "Ancestry", "Age", "Sex"),
    assoc_list = assoc_list,
    col_rows = list(
        'APOE_carrier' = APOE_carrier_colors,
        'Ancestry' = ancestry_colors,
        'Sex' = sex_colors
    )
)
dev.off()

#### MOFA gene weights ####

rd <- rowData(spe) |>
    as.data.frame() |>
    select(feature = gene_id, gene_name)

factor_names <- paste0("Factor", seq(model_opts$num_factors))
names(factor_names) <- factor_names

gene_weights <- map(factor_names, ~MOFAcellulaR::get_geneweights(model = model, factor = .x) |> left_join(rd, by = join_by(feature)))
map(gene_weights, dim)

write_rds(gene_weights, file = here(data_dir, sprintf('MOFA_gene_weights_%s.pdf', opt$datatype)))

# slurmjobs::job_single('01_MOFA_broad', create_shell = TRUE, memory = '10G', command = "Rscript 04_DEG_boxplots.R --datatype sn_broad")
# slurmjobs::job_single('01_MOFA_fine', create_shell = TRUE, memory = '10G', command = "Rscript 04_DEG_boxplots.R --datatype sn_fine")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()