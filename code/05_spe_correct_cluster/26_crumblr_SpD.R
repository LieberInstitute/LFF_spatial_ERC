## Louise Huuki-Myers, April 2025
## Run differential proportion analysis with Crumblr

library("SingleCellExperiment")
library("tidyverse")
library("crumblr")
library("variancePartition")
library("here")
library("sessioninfo")
library("dendextend")
library("jaffelab")

## Prep directories
plot_dir <- here("plots", "05_spe_correct_cluster", "26_crumblr_SpD")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "05_spe_correct_cluster", "26_crumblr_SpD")
if(!dir.exists(data_dir)) dir.create(data_dir)

## load colors
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
load(here("processed-data", "05_spe_correct_cluster", "19_SpD_update_spe", "SpD_proportions.Rdata"), verbose = TRUE)
# SpD_proportions

SpD_proportions <- SpD_proportions |> filter(sample_id != "Br1289")

SpD_df <- SpD_proportions |> 
    mutate(SpD_syn = gsub("~", "_", SpD)) |>
    ungroup() |>
    select(SpD, SpD_syn) |>
    unique()

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- SpD_proportions |>
    mutate(error = sqrt((prop * (1-prop))/sum(n)),
           anno = paste0(n,"/", sum(n))) |>
    select(sample_id, SpD, error, anno)

## create cell count matrix
cell_counts <- SpD_proportions |> 
    mutate(SpD_syn = gsub("~", "_", SpD)) |>
    ungroup() |>
    select(sample_id, SpD_syn, n) |>
    pivot_wider(names_from = SpD_syn, values_from = n)  |>
    column_to_rownames("sample_id")

## replace NA w/ 0
any(is.na(cell_counts))
cell_counts[is.na(cell_counts)] <- 0

## create info matrix
erc_info <- read.csv(here("processed-data", "02_build_spe", "sample_info.csv")) |> 
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE, Visium_slide, round) |>
    mutate(APOE_carrier = ifelse(grepl("E2", APOE), "E2+", "E4+")) |>
    filter(sample_id != "Br1289") |>
    column_to_rownames("sample_id")

setequal(rownames(erc_info), rownames(cell_counts))

erc_info <- erc_info[rownames(cell_counts), ]

identical(rownames(erc_info), rownames(cell_counts))

#### Apply crumblr transformation ####
cobj <- crumblr(cell_counts)

# E:numeric matrix of CLR transformed counts

clr_prop_long <- cobj$E |>
    as.data.frame() |>
    rownames_to_column("SpD_syn") |>
    pivot_longer(!SpD_syn, names_to = "sample_id", values_to = "CLR") |>
    left_join(erc_info |> rownames_to_column("sample_id")) |>
    left_join(SpD_df) |>
    left_join(error)

# Partition variance into components for sample
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | Visium_slide)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

ggsave(fig.vp + 
           theme_bw() + 
           theme(legend.position = "bottom"),
       filename = here(plot_dir, "crumblr_SpD_vp.png"), width = 4.5, height = 5)

## variable correlation
C <- canCorPairs(form, erc_info)
# APOE and APOE_carrier
# Visium_slide and round
pdf(here(plot_dir, "variable_CorrMatrix_SpD.pdf"), height = 5, width = 5)
plotCorrMatrix(C)
dev.off()

#### PCA ####
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
#   shape by Stimulated vs unstimulated
cobj_pca <- ggplot(df_pca, aes(PC1, PC2, color = Visium_slide, shape = APOE)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

ggsave(cobj_pca, filename = here(plot_dir, "crumblr_pca_SpD.png"))

#### CLR plots ####

clr_boxplot_APOE <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~SpD) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE, filename = here(plot_dir, "clr_boxplot_APOE.png"))

clr_boxplot_APOE_carrier <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~SpD, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier, filename = here(plot_dir, "clr_boxplot_APOE_carrier.png"), height = 3, width = 10)

clr_boxplot_Visium_slide <- clr_prop_long |>
    ggplot(aes(x = Visium_slide, y = CLR, fill = Visium_slide)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~SpD) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Visium_slide, filename = here(plot_dir, "clr_boxplot_Visium_slide.png"))


#### Hierarchical clustering ####
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "25_SpD_pseudobulk", "spe_pseudobulk_only-SpD_syn.rds"))
dim(spe_pb)

dist.clusCollapsed <- dist(t(logcounts(spe_pb)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)
# plot(dend, main = "hierarchical cluster dend", horiz = TRUE)

# Use variancePartition workflow to analyze each SpD

#### fit APOE_carrier ####
fit <- dream(cobj, ~ APOE_carrier + Anc_Afr + Age + Visium_slide, erc_info)
fit <- eBayes(fit)

# Test APOE carrier
(diff_prop_APOE_carrier <- topTable(fit, coef = "APOE_carrierE4+", number = Inf))
#                     logFC     AveExpr           t    P.Value adj.P.Val         B
# WM_Sp09D06    -0.49369667 -0.06336945 -1.88330104 0.07319348 0.6587413 -4.276222
# L2.3_Sp09D01   0.29782763 -0.21215444  1.36020054 0.18779244 0.7776110 -4.487527
# Inhib_Sp09D09  0.11282175 -0.51510636  0.97529404 0.34020990 0.7776110 -4.627923
# L6_Sp09D04    -0.10189724  0.55858775 -0.87316693 0.39217475 0.7776110 -4.645685
# L5_Sp09D03    -0.13401893 -0.16236857 -0.68279468 0.50199373 0.7776110 -4.759007
# WM.uf_Sp09D07  0.15135687 -0.06844375  0.49499131 0.62560326 0.7776110 -4.733619
# LD_Sp09D02     0.04024377  0.85387084  0.43449588 0.66823749 0.7776110 -4.723201
# L1_Sp09D05    -0.09734721  0.28070405 -0.40258468 0.69120979 0.7776110 -4.702966
# Vasc_Sp09D08  -0.01867111 -0.67172007 -0.07730314 0.93909390 0.9390939 -4.762624

write.csv(diff_prop_APOE_carrier, file = here(data_dir, "SpD_diff_prop_APOE_carrier.csv"))

#### cleanY on cobj  APOE_carrier ####

mod <- model.matrix( ~ APOE_carrier + Sex + Age + Anc_Afr + Visium_slide , erc_info)
cleanY_cobj <- cleaningY(cobj$E , mod, P=2)
head(cleanY_cobj)

clean_clr_prop_long <- cleanY_cobj |>
    as.data.frame() |>
    rownames_to_column("SpD_syn") |>
    pivot_longer(!SpD_syn, names_to = "sample_id", values_to = "CLR_cleanY")

clr_prop_long <- clr_prop_long |> left_join(clean_clr_prop_long)

## clr plot of top results
clr_cleanY_boxplot_APOE_carrier  <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~SpD, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_cleanY_boxplot_APOE_carrier, filename = here(plot_dir, "clr_cleanY_boxplot_APOE_carrier.png"), height = 4, width = 12)


# Perform multivariate test across the hierarchy - APOE_carrier
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

res_tb <- res |> as_tibble()
res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "SpD_diff_prop_tree_test_APOE_carrier.csv"))

# Plot hierarchy and testing results
tree_fdr_plot <- plotTreeTest(res)
ggsave(tree_fdr_plot, filename = here(plot_dir, "crumblr_SpD_APOE_carrier_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
tree_beta_plot <- plotTreeTestBeta(res, low = "#398A84", high = "#D46B43") +
    theme(legend.position = "bottom", 
          legend.box = "vertical"
    )
ggsave(tree_beta_plot, filename = here(plot_dir, "crumblr_SpD_APOE_carrier_Tree_beta.png"), width = 4, height = 5)

# forest plots
forest_plot <- plotForest(res, hide = FALSE, low = "#398A84", high = "#D46B43") 
ggsave(forest_plot, filename = here(plot_dir, "crumblr_SpD_APOE_carrier_forest.png"), width = 4, height = 5)


# combined_fig <- fig.vp +
#     theme(legend.position = "left") |
#     plotTreeTestBeta(res) +
#     theme(legend.position = "bottom", legend.box = "vertical") |
#     plotForest(res, hide = FALSE) 
# 
# ggsave(combined_fig, filename = here(plot_dir, "crumblr_SpD_combined_APOE_carrier.png"))

# Perform multivariate test across the hierarchy - visium slide
(slides <- colnames(fit$design)[grep("Visium_slide", colnames(fit$design))])

topTable(fit, coef = slides, number = Inf)
#               Visium_slideV13Y24.343 Visium_slideV13Y24.344     AveExpr         F    P.Value adj.P.Val
# L3_Sp09D02               -0.09108499             -0.1492597  0.85933332 3.2954274 0.01356117 0.1220505
# WM_Sp09D06                0.28476468             -0.2938779 -0.06954456 1.3910685 0.25485856 0.5457378
# WM.uf_Sp09D07             0.33007275              0.2088800 -0.08853475 1.3733007 0.26193527 0.5457378
# L6_Sp09D04               -0.04727934              0.0838746  0.57643431 1.2683854 0.30746514 0.5457378
# L2.3_Sp09D01              0.05055732             -0.5125162 -0.19707896 1.1430306 0.37079949 0.5457378
# Vasc_Sp09D08              0.59171416              0.6716892 -0.70635115 1.0920655 0.39946856 0.5457378
# L5_Sp09D03                0.15539167             -0.2040875 -0.14699867 1.0362420 0.43283429 0.5457378
# L1_Sp09D05                0.48705349              0.4129016  0.27080407 0.9548491 0.48510030 0.5457378
# LD_Sp09D09               -0.10479198             -0.1859432 -0.49806362 0.8272432 0.57497078 0.5749708

# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "FE")
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "FE.empirical")
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "RE2C")
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "tstat")
## only these work for muliple covar
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "sidak")
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "fisher")

# Plot hierarchy and testing results
plotTreeTest(res)
plotTreeTestBeta(res) 

plotForest(res)
           
combined_fig <- fig.vp +
    theme(legend.position = "left") |
    plotTreeTestBeta(res) +
    theme(legend.position = "bottom", legend.box = "vertical") |
    plotForest(res, hide = FALSE) 

ggsave(combined_fig, filename = here(plot_dir, "crumblr_SpD_combined_APOE_carrier.png"))


#### fit APOE_carrier ####
fit <- dream(cobj, ~ APOE + Age + Visium_slide, erc_info)
fit <- eBayes(fit)

# Test APOE carrier
(geno <- colnames(fit$design)[grep("APOE", colnames(fit$design))])

topTable(fit, coef = geno, number = Inf)
#                 APOEE2.E3   APOEE3.E4   APOEE4.E4     AveExpr         F    P.Value  adj.P.Val
# WM.uf_Sp09D07 -1.45822296 -0.71161195 -0.57578954 -0.08853475 4.7981628 0.01019429 0.07254207
# L2.3_Sp09D01   0.97780063  1.05490717  0.91721907 -0.19707896 4.2715653 0.01612046 0.07254207
# L1_Sp09D05    -0.98712303 -0.68638232 -0.69236950  0.27080407 3.5695900 0.03053114 0.09159341
# WM_Sp09D06    -0.72657490 -0.91152411 -0.83819807 -0.06954456 2.6446709 0.07451405 0.14996941
# LD_Sp09D09     0.43800262  0.49694333  0.40256131 -0.49806362 2.5329896 0.08331634 0.14996941
# L5_Sp09D03     0.62567182  0.37905216  0.27819314 -0.14699867 1.6860349 0.19918238 0.29877357
# L3_Sp09D02     0.22188381  0.24068605  0.11585516  0.85933332 1.0709642 0.38178055 0.49086070
# Vasc_Sp09D08  -0.13047354  0.03866813 -0.34979559 -0.70635115 0.6093153 0.61607566 0.69308512
# L6_Sp09D04     0.03480808 -0.11499476  0.05244457  0.57643431 0.3705864 0.77498382 0.77498382


topTable(fit, coef = "APOEE2/E3", number = Inf)
#                     logFC     AveExpr          t     P.Value  adj.P.Val          B
# WM.uf_Sp09D07 -1.45822296 -0.08853475 -3.7192489 0.001198532 0.01078679 -0.9108791
# L1_Sp09D05    -0.98712303  0.27080407 -3.1978520 0.004165307 0.01435555 -1.9918766
# L2.3_Sp09D01   0.97780063 -0.19707896  3.1387595 0.004785185 0.01435555 -2.3535928
# LD_Sp09D09     0.43800262 -0.49806362  2.3705441 0.026979388 0.06070362 -3.7208538
# L5_Sp09D03     0.62567182 -0.14699867  2.2111866 0.037757567 0.06785159 -4.2650056
# WM_Sp09D06    -0.72657490 -0.06954456 -2.1235449 0.045234391 0.06785159 -4.3347318
# L3_Sp09D02     0.22188381  0.85933332  1.5814708 0.128090230 0.16468744 -5.0093801
# Vasc_Sp09D08  -0.13047354 -0.70635115 -0.3812269 0.706704178 0.79504220 -6.3409796
# L6_Sp09D04     0.03480808  0.57643431  0.1756788 0.862158770 0.86215877 -6.1948897

# Perform multivariate test across the hierarchy - APOE geno
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOEE2/E3")

# Plot hierarchy and testing results
plotTreeTest(res)
plotTreeTestBeta(res) 
plotForest(res)


# slurmjobs::job_single('26_crumblr_SpD', create_shell = TRUE, memory = '25G', command = "Rscript 26_crumblr_SpD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
