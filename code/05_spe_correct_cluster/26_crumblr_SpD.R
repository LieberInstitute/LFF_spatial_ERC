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

SpD_levels <- names(SpD_colors)

#### Load data ####
load(here("processed-data", "05_spe_correct_cluster", "19_SpD_update_spe", "SpD_proportions.Rdata"), verbose = TRUE)
# SpD_proportions

SpD_proportions <- SpD_proportions |> filter(sample_id != "Br1289")

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- SpD_proportions |>
    mutate(error = sqrt((prop * (1-prop))/sum(n)),
           anno = paste0(n,"/", sum(n))) |>
    select(sample_id, vSpD, error, anno)

## create cell count matrix
cell_counts <- SpD_proportions |> 
    ungroup() |>
    select(sample_id, vSpD, n) |>
    pivot_wider(names_from = vSpD, values_from = n)  |>
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
    rownames_to_column("vSpD") |>
    pivot_longer(!vSpD, names_to = "sample_id", values_to = "CLR") |>
    left_join(erc_info |> rownames_to_column("sample_id")) |>
    left_join(error)

# Partition variance into components for sample
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | Visium_slide)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
# fig.vp

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
    # scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

ggsave(cobj_pca, filename = here(plot_dir, "crumblr_pca_SpD.png"))

#### CLR plots ####

clr_boxplot_APOE <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~vSpD) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE, filename = here(plot_dir, "clr_boxplot_APOE.png"))

clr_boxplot_APOE_carrier <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~vSpD, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier, filename = here(plot_dir, "clr_boxplot_APOE_carrier.png"), height = 3, width = 10)

clr_boxplot_Visium_slide <- clr_prop_long |>
    ggplot(aes(x = Visium_slide, y = CLR, fill = Visium_slide)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~vSpD) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Visium_slide, filename = here(plot_dir, "clr_boxplot_Visium_slide.png"))


#### Hierarchical clustering ####
spe_pb <- readRDS(here("processed-data", "05_spe_correct_cluster", "25_SpD_pseudobulk", "spe_pseudobulk_only-vSpD.rds"))
dim(spe_pb)

dist.clusCollapsed <- dist(t(logcounts(spe_pb)))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)

pdf(here(plot_dir, "hierarchical_cluster_dendrogram_vSpD.pdf"))
plot(dend, main = "hierarchical cluster dend", horiz = TRUE)
dev.off()

# Use variancePartition workflow to analyze each SpD

#### fit APOE_carrier ####
fit <- dream(cobj, ~ APOE_carrier + Anc_Afr + Age + Visium_slide, erc_info)
fit <- eBayes(fit)

# Test APOE carrier
(diff_prop_APOE_carrier <- topTable(fit, coef = "APOE_carrierE4+", number = Inf))
#              logFC    AveExpr           t    P.Value adj.P.Val         B
# vWMd   -0.69198449 -0.7925698 -2.30943789 0.03000760 0.3300836 -3.524001
# vL2     0.48733930  0.2590005  1.79309029 0.08580174 0.3959066 -4.118591
# vL5b    0.47879560  0.2865521  1.67086840 0.10797452 0.3959066 -4.247660
# vInhib  0.20547960 -1.2511731  1.32306758 0.19850529 0.5458896 -4.527167
# vWMpv  -0.32388097 -0.1606711 -1.05357860 0.30275696 0.6280980 -4.757622
# vL6    -0.17906856  0.5290294 -0.96856016 0.34259889 0.6280980 -4.813384
# vL5a   -0.17744437  0.1263453 -0.74014624 0.46652027 0.6674069 -4.978527
# vL1    -0.21874928  0.6487558 -0.70882259 0.48538681 0.6674069 -4.834797
# vLD     0.10152290  0.7206322  0.58551407 0.56376658 0.6779896 -4.931778
# vWMuf   0.17037805  0.1973571  0.50773999 0.61635418 0.6779896 -4.993951
# vVasc  -0.01696857 -0.5632584 -0.05991562 0.95272795 0.9527280 -5.026387

write.csv(diff_prop_APOE_carrier, file = here(data_dir, "SpD_diff_prop_APOE_carrier.csv"))

#### cleanY on cobj  APOE_carrier ####

mod <- model.matrix( ~ APOE_carrier + Sex + Age + Anc_Afr + Visium_slide , erc_info)
cleanY_cobj <- cleaningY(cobj$E , mod, P=2)
head(cleanY_cobj)

clean_clr_prop_long <- cleanY_cobj |>
    as.data.frame() |>
    rownames_to_column("vSpD") |>
    pivot_longer(!vSpD, names_to = "sample_id", values_to = "CLR_cleanY")

clr_prop_long <- clr_prop_long |> 
    left_join(clean_clr_prop_long) |>
    mutate(vSpD = factor(vSpD, levels = SpD_levels))

## clr plot of top results
clr_cleanY_boxplot_APOE_carrier  <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~vSpD, scales = "free_y", nrow = 1) +
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
# Visium_slideV13Y24.342 Visium_slideV13Y24.343 Visium_slideV13Y24.344    AveExpr         F    P.Value adj.P.Val
# vWMd              -0.80351260             0.57154893           -0.898452477 -0.7925698 2.7067623 0.03274812 0.3602293
# vLD                0.09468673            -0.33180266           -0.140153068  0.7206322 1.7128547 0.15432382 0.6206414
# vL2                0.56641285             0.22150738           -0.444922193  0.2590005 1.3405776 0.27571861 0.6206414
# vL1                0.37281046             0.75339193            0.657143517  0.6487558 1.2944525 0.29580661 0.6206414
# vL6                0.27671656             0.09665595            0.255609798  0.5290294 1.2376391 0.32232400 0.6206414
# vWMuf             -0.47356866             0.24001734            0.136346822  0.1973571 1.2048821 0.33853165 0.6206414
# vL5a               0.62905803             0.84072543            0.326752297  0.1263453 0.9079456 0.51718765 0.7418922
# vL5b               0.23809129            -0.24496856           -0.289252154  0.2865521 0.8669803 0.54612801 0.7418922
# vVasc              0.68481767             0.66452535            0.850362758 -0.5632584 0.7844065 0.60700271 0.7418922
# vWMpv             -0.42367786             0.52231617            0.009814914 -0.1606711 0.5770769 0.76733878 0.7877749
# vInhib             0.12542752            -0.09287964           -0.068536604 -1.2511731 0.5502414 0.78777489 0.7877749

# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "FE")
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "FE.empirical")
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "RE2C")
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "tstat")
## only these work for muliple covar
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "sidak")
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = slides, method = "fisher")

# Plot hierarchy and testing results
# plotTreeTest(res)
# plotTreeTestBeta(res) 
# 
# plotForest(res)
#            
# combined_fig <- fig.vp +
#     theme(legend.position = "left") |
#     plotTreeTestBeta(res) +
#     theme(legend.position = "bottom", legend.box = "vertical") |
#     plotForest(res, hide = FALSE) 
# 
# ggsave(combined_fig, filename = here(plot_dir, "crumblr_SpD_combined_APOE_carrier.png"))
# 

#### fit APOE_carrier ####
fit <- dream(cobj, ~ APOE + Age + Visium_slide, erc_info)
fit <- eBayes(fit)

# Test APOE carrier
(geno <- colnames(fit$design)[grep("APOE", colnames(fit$design))])

topTable(fit, coef = geno, number = Inf)
#         APOEE2.E3    APOEE3.E4   APOEE4.E4    AveExpr         F    P.Value adj.P.Val
# vWMuf  -1.4792312 -0.678682931 -0.43884724  0.1973571 4.4893064 0.01282457 0.1272176
# vL1    -1.2232885 -1.011312475 -0.74114774  0.6487558 3.4198735 0.03431942 0.1272176
# vL2     1.0675820  1.323825639  1.02873329  0.2590005 3.4084638 0.03469570 0.1272176
# vWMd   -0.6416058 -1.024912066 -1.02029635 -0.7925698 2.3714389 0.09689224 0.2551571
# vL5a    0.7958566  0.302081456  0.47584396  0.1263453 2.1513275 0.12155455 0.2551571
# vWMpv  -0.9179539 -0.849297602 -0.77556945 -0.1606711 2.0211745 0.13917659 0.2551571
# vL5b    0.7563362  1.087187775  0.68640148  0.2865521 1.7372074 0.18755243 0.2947252
# vLD     0.1136103  0.421578891 -0.03805594  0.7206322 1.4608008 0.25149517 0.3458059
# vInhib  0.3620315  0.464401527  0.36140096 -1.2511731 1.1559260 0.34805307 0.4253982
# vL6    -0.1240307 -0.351651950 -0.20753701  0.5290294 0.5464017 0.65556579 0.7211224
# vVasc  -0.1677176  0.003933727 -0.33546314 -0.5632584 0.3436182 0.79402577 0.7940258


topTable(fit, coef = "APOEE2/E3", number = Inf)
#             logFC    AveExpr          t     P.Value  adj.P.Val         B
# vWMuf  -1.4792312  0.1973571 -3.4975452 0.001954218 0.02149640 -1.332302
# vL1    -1.2232885  0.6487558 -3.1430477 0.004582582 0.02520420 -2.040101
# vL2     1.0675820  0.2590005  2.6258585 0.015157620 0.05557794 -3.258967
# vL5a    0.7958566  0.1263453  2.3184468 0.029734716 0.07059277 -3.871608
# vWMpv  -0.9179539 -0.1606711 -2.2826970 0.032087625 0.07059277 -3.963263
# vL5b    0.7563362  0.2865521  1.6991708 0.102868408 0.18859208 -4.933679
# vWMd   -0.6416058 -0.7925698 -1.5233722 0.141392049 0.20699618 -5.250391
# vInhib  0.3620315 -1.2511731  1.4875615 0.150542675 0.20699618 -5.042012
# vLD     0.1136103  0.7206322  0.4664093 0.645340335 0.67946447 -6.016474
# vL6    -0.1240307  0.5290294 -0.4465423 0.659410553 0.67946447 -6.175231
# vVasc  -0.1677176 -0.5632584 -0.4185361 0.679464469 0.67946447 -6.085853

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
