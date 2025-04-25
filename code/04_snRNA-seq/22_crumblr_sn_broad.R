## Louise Huuki-Myers, April 2025
## Run differential proportion analysis with Crumblr

library("SingleCellExperiment")
library("tidyverse")
library("crumblr")
library("variancePartition")
library("dendextend")
library("here")
library("sessioninfo")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "22_crumblr_sn_broad")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "22_crumblr_sn_broad")
if(!dir.exists(data_dir)) dir.create(data_dir)

## load colors
load(here("processed-data", "04_snRNA-seq", "cell_type_colors.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
# cell_type_proportions
load(here("processed-data", "04_snRNA-seq", "15_cluster_update_sce","cell_type_proportions.Rdata"), verbose = TRUE)

## calc broad proportions 
cell_type_proportions <- cell_type_proportions |>
    separate(cell_type_anno, into = c("cell_type_broad"), sep ="\\.", extra = "drop") |>
    mutate(cell_type_broad = factor(cell_type_broad, levels = c("Astro","Endo","Micro","Oligo","OPC","Excit","Inhib"))) |>
    group_by(sample_id, APOE, Sex, Age, Ancestry, Anc_Afr, exp_round, seq_round, cell_type_broad) |>
    summarise(n = sum(n), prop = sum(prop)) |>
    group_by(sample_id)

# cell_type_proportions_broad |> summarise(sum(prop))

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- cell_type_proportions |> #already grouped by sample
    mutate(error = sqrt((prop * (1-prop))/sum(n)),
           anno = paste0(n,"/", sum(n))) |>
    select(sample_id, cell_type_broad, prop, error, anno)

## create cell count matrix
cell_counts <- cell_type_proportions |> 
    ungroup() |>
    select(sample_id, cell_type_broad, n) |>
    pivot_wider(names_from = cell_type_broad, values_from = n)  |>
    column_to_rownames("sample_id")

## replace NA w/ 0
cell_counts[is.na(cell_counts)] <- 0

## create info matrix

erc_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv")) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, exp_round, seq_round)|>
    column_to_rownames("sample_id")

setequal(rownames(erc_info), rownames(cell_counts))

erc_info <- erc_info[rownames(cell_counts), ]
identical(rownames(erc_info), rownames(cell_counts))

#### Apply crumblr transformation ####
cobj <- crumblr(cell_counts)

# E:numeric matrix of CLR transformed counts
clr_prop_long <- cobj$E |>
    as.data.frame() |>
    rownames_to_column("cell_type_broad") |>
    pivot_longer(!cell_type_broad, names_to = "sample_id", values_to = "CLR") |>
    left_join(erc_info |> rownames_to_column("sample_id")) |>
    left_join(error) |>
    replace_na(list(prop = 0, error = 0))

# Partition variance into components for sample
form <- ~ (1 | APOE) + (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | exp_round)+ (1 | seq_round)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

ggsave(fig.vp, filename = here(plot_dir, "crumblr_cell_type_broad_vp.png"))

#### variable correlation ####
C <- canCorPairs(form, erc_info)
# APOE and APOE_carrier
pdf("variable_CorrMatrix_sn_broad.pdf", height = 5, width = 5)
plotCorrMatrix(C)
dev.off()

#### PCA ####
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
#   shape by Stimulated vs unstimulated
ggplot(df_pca, aes(PC1, PC2, color = exp_round, shape = APOE)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

#### CLR plots ####

clr_v_prop <- clr_prop_long |>
    ggplot(aes(prop, CLR, color = error)) +
    geom_point()

clr_boxplot_APOE <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot() +
    facet_wrap(~cell_type_broad) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE, filename = "clr_boxplot_APOE_cell_type_broad.png")

clr_boxplot_APOE_jitter <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_broad) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_jitter, filename = "clr_boxplot_APOE_jitter_cell_type_broad.png")

clr_boxplot_APOE_carrier <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_broad) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_carrier, filename = "clr_boxplot_APOE_carrier_cell_type_broad.png")

clr_boxplot_exp_round <- clr_prop_long |>
    ggplot(aes(x = exp_round, y = CLR, fill = exp_round)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_broad) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_exp_round, filename = "clr_boxplot_exp_round.png")

#### CLR + prop plots ####

clr_prop_long2 <- clr_prop_long |>
    pivot_longer(cols = c("CLR", "prop"), names_to = "Metric")


clr_prop_long2 |>
    ggplot(aes(x = APOE, y = value, fill = APOE)) +
    geom_boxplot() +
    facet_grid(Metric~cell_type_broad, scales = "free") +
    # facet_wrap(cell_type_broad~Metric, scales = "free", ncol = 2) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()


#### Hierarchical clustering ####
sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_pseudobulk_only-cell_type_broad.rds"))
dim(sce_pb)

dist.clusCollapsed <- dist(t(logcounts(sce_pb)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)
# plot(dend, main = "hierarchical cluster dend", horiz = TRUE)

#### fit APOE_carrier ####
# Use variancePartition workflow to analyze each cell type
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance

fit <- dream(cobj, ~ APOE_carrier + Sex + Age + Anc_Afr + exp_round , erc_info)
fit <- eBayes(fit)

# Extract results for each cell type
(diff_prop_APOE_carrier <- topTable(fit, coef = "APOE_carrierE4+", number = Inf))
#             logFC     AveExpr           t     P.Value   adj.P.Val         B
# Inhib  0.405945782 -0.09286424  3.61177610 0.001381888 0.009673219 -1.056731
# Oligo -0.477374911  1.36138934 -2.69783115 0.012513726 0.043798042 -3.136641
# Excit  0.230026614  0.06625368  1.54241900 0.135940684 0.317194928 -5.287135
# OPC    0.053545961  0.26262646  0.37642310 0.709880619 0.989634045 -6.166996
# Astro -0.040679217  0.60484917 -0.34592890 0.732383741 0.989634045 -6.270353
# Endo  -0.004801482 -2.24991564 -0.02133842 0.983150564 0.989634045 -6.400113
# Micro -0.002069250  0.04766122 -0.01312698 0.989634045 0.989634045 -6.497953

write.csv(diff_prop_APOE_carrier, file = here(data_dir, "diff_prop_cell_type_broad_APOE_carrier.csv"))

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

# Plot hierarchy and testing results
plotTreeTest(res)

# Plot hierarchy and regression coefficients
plotTreeTestBeta(res)
plotForest(res, hide = FALSE) 

combined_fig <- fig.vp +
    theme(legend.position = "left") |
    plotTreeTestBeta(res) +
    theme(legend.position = "bottom", legend.box = "vertical") |
    plotForest(res, hide = FALSE) 

ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined_APOE_carrier.png"), width = 10)

## clr plot of top results
clr_boxplot_APOE_carrier_Excit <- clr_prop_long |>
    filter(cell_type_broad == "Excit") |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_broad) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_carrier_Excit, filename = here(plot_dir, "clr_boxplot_APOE_carrier_Excit.png"))

prop_boxplot_APOE_carrier_Excit <- clr_prop_long |>
    filter(cell_type_broad == "Excit") |>
    ggplot(aes(x = APOE_carrier, y = prop, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_broad) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(prop_boxplot_APOE_carrier_Excit, filename = here(plot_dir, "prop_boxplot_APOE_carrier_Excit.png"))

#### fit exp_round ####
# Extract results for each cell type
topTable(fit, coef = "exp_roundround2", number = Inf)

#                        logFC    AveExpr           t   P.Value adj.P.Val         B
# Astro.1          -0.61245101  2.1654993 -1.55953049 0.1328783 0.9693659 -4.573971
# Inhib.Lamp5_Lhx6  0.61170837 -0.4007870  1.48596918 0.1512167 0.9693659 -4.574353
# Excit.L2         -0.95274566 -0.1314255 -1.25339076 0.2229793 0.9693659 -4.588706
# Inhib.Pvalb       0.63602081 -0.4311626  1.24669577 0.2253792 0.9693659 -4.582203

topTable(fit, coef = "SexM", number = Inf)

# write.csv(diff_prop_APOE_carrier, file = here(data_dir, "diff_prop_APOE_carrier.csv"))


#### fit APOE ####

fit2 <- dream(cobj, ~ APOE + Age + Anc_Afr + exp_round , erc_info)
fit2 <- eBayes(fit2)

# Test APOE carrier
(geno <- colnames(fit2$design)[grep("APOE", colnames(fit2$design))])
(diff_prop_APOE <- topTable(fit2, coef = geno, number = Inf))
#                   APOEE2.E3    APOEE3.E4    APOEE4.E4    AveExpr         F    P.Value adj.P.Val
# Excit.L6b        -0.3774346 -1.642397391 -1.017496429 -1.0645145 4.7032782 0.01140097 0.2374917
# Oligo.2          -0.3884539 -1.142631475 -0.882662486  2.2535239 3.7707442 0.02588605 0.2374917
# Excit.L2          0.1743530  0.984704658  0.653632584 -0.1314255 3.2310558 0.04272876 0.2374917

write.csv(diff_prop_APOE, file = here(data_dir, "diff_prop_APOE.csv"))


# Perform multivariate test across the hierarchy
topTable(fit, coef = "APOEE4/E4", number = Inf)

res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOEE4/E4")

# Plot hierarchy and testing results
plotTreeTest(res)

# Plot hierarchy and regression coefficients
plotTreeTestBeta(res)
plotForest(res, hide = FALSE) 

combined_fig <- fig.vp +
    theme(legend.position = "left") |
    plotTreeTestBeta(res) +
    theme(legend.position = "bottom", legend.box = "vertical") |
    plotForest(res, hide = FALSE) 

ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined_APOEE4.E4.png"))

# slurmjobs::job_single('22_crumblr_sn_broad', create_shell = TRUE, memory = '25G', command = "Rscript 19_sn_heatmaps.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
