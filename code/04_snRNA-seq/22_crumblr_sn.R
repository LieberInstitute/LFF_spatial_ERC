## Louise Huuki-Myers, April 2025
## Run differential proportion analysis with Crumblr

# BiocManager::install('DiseaseNeurogenomics/crumblr')

library("SingleCellExperiment")
library("tidyverse")
library("crumblr")
library("variancePartition")
library("dendextend")
library("here")
library("sessioninfo")

## Prep directories
plot_dir <- here("plots", "04_snRNA-seq", "22_crumblr_sn")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "04_snRNA-seq", "22_crumblr_sn")
if(!dir.exists(data_dir)) dir.create(data_dir)

## load colors
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
# cell_type_proportions
load(here("processed-data", "04_snRNA-seq", "28_subcluster_update_sce","cell_type_proportions.Rdata"), verbose = TRUE)

cell_type_proportions <- cell_type_proportions |> filter(sample_id != "Br1289")

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- cell_type_proportions |> #already grouped by sample
    mutate(error = sqrt((prop * (1-prop))/sum(n)),
           anno = paste0(n,"/", sum(n))) |>
    select(sample_id, cell_type_anno, prop, error, anno)

## create cell count matrix
cell_counts <- cell_type_proportions |> 
    ungroup() |>
    select(sample_id, cell_type_anno, n) |>
    pivot_wider(names_from = cell_type_anno, values_from = n)  |>
    column_to_rownames("sample_id")

## replace NA w/ 0
cell_counts[is.na(cell_counts)] <- 0

## create info matrix

erc_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv")) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE, exp_round, seq_round) |> 
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
    rownames_to_column("cell_type_anno") |>
    pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR") |>
    left_join(erc_info |> rownames_to_column("sample_id")) |>
    left_join(error) |>
    replace_na(list(prop = 0, error = 0))

# Partition variance into components for sample
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | exp_round)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

ggsave(fig.vp + 
           theme_bw() + 
           theme(legend.position = "bottom"),
       filename = here(plot_dir, "crumblr_cell_type_vp.png"), width = 4.5, height = 7)

#### variable correlation ####
C <- canCorPairs(form, erc_info)
# APOE and APOE_carrier
pdf(here(plot_dir, "variable_CorrMatrix_sn.pdf"), height = 5, width = 5)
plotCorrMatrix(C)
dev.off()

#### PCA ####
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
#   shape by Stimulated vs unstimulated
cobj_pca <- ggplot(df_pca, aes(PC1, PC2, color = exp_round, shape = APOE)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

ggsave(cobj_pca, filename = here(plot_dir, "crumblr_pca_sn.png"))

#### CLR plots ####

clr_prop_long |> filter(cell_type_anno == "Oligo.3") |> arrange(-prop)

clr_v_prop <- clr_prop_long |>
    ggplot(aes(prop, CLR, color = error)) +
    geom_point()

clr_boxplot_APOE <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot() +
    facet_wrap(~cell_type_anno) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE, filename = here(plot_dir, "clr_boxplot_APOE.png"))

clr_boxplot_APOE_jitter <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_anno) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_jitter, filename = here(plot_dir, "clr_boxplot_APOE_jitter.png"))

clr_boxplot_APOE_carrier <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_anno) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_carrier, filename = here(plot_dir, "clr_boxplot_APOE_carrier.png"))

clr_boxplot_exp_round <- clr_prop_long |>
    ggplot(aes(x = exp_round, y = CLR, fill = exp_round)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_anno) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_exp_round, filename = here(plot_dir, "clr_boxplot_exp_round.png"))


#### Hierarchical clustering ####
# sce_pb <- readRDS(here("processed-data", "04_snRNA-seq", "20_sn_pseudobulk","sce_pseudobulk_only-cell_type_anno.rds"))
# dim(sce_pb)
# 
# dist.clusCollapsed <- dist(t(logcounts(sce_pb)))
# tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")
# 
# dend <- as.dendrogram(tree.clusCollapsed, hang = 0.2)
# plot(dend, main = "hierarchical cluster dend", horiz = TRUE)

load(here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster","sn_subcluster_hierarchical_cluster.Rdata"), verbose = TRUE)
# dend
# tree.clusCollapsed
# dist.clusCollapsed


#### Run DREAM + eBayes on cobj ####
fit <- dream(cobj, ~ APOE_carrier + Sex + Age + Anc_Afr + exp_round , erc_info) 
fit <- eBayes(fit = fit)

#### fit APOE_carrier ####
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance

# Extract results for each cell type
(diff_prop_APOE_carrier <- topTable(fit, coef = "APOE_carrierE4+", number = Inf))

#                         logFC     AveExpr           t     P.Value  adj.P.Val         B
# Oligo.1          -0.838001493  2.12294822 -3.29041883 0.003391148 0.07677455 -1.714310
# Excit.L2          0.878852848  0.06115206  3.21587912 0.004040766 0.07677455 -1.981697
# Oligo.2          -0.903751212  1.90335046 -3.02138982 0.006355225 0.07762043 -2.274637
# Excit.L6b        -1.093021726 -0.88692098 -2.91208013 0.008170572 0.07762043 -2.577251
# Inhib.Pvalb       0.697314283 -0.23355939  2.77258192 0.011215030 0.08523423 -2.932870
# Inhib.Lamp5_Lhx6  0.546190860 -0.18985351  2.63405746 0.015284600 0.09680246 -3.112341
# Astro.4          -0.707230347 -0.37611093 -2.32148552 0.030091782 0.15319537 -3.821864
# OPC.5             1.088450759 -1.80267586  2.28861479 0.032251657 0.15319537 -3.891058

write.csv(diff_prop_APOE_carrier, file = here(data_dir, "diff_prop_APOE_carrier.csv"))

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "sn_diff_prop_tree_test_APOE_carrier.csv"))

# Plot hierarchy and testing results
tree_fdr_plot <- plotTreeTest(res)
ggsave(tree_fdr_plot, filename = here(plot_dir, "crumblr_cell_type_APOE_carrier_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
# TODO fdr.cutoff = 0.1
tree_beta_plot <- plotTreeTestBeta(res, low = "#398A84", high = "#D46B43") +
    theme(legend.position = "left" 
          # legend.box = "vertical"
    )

ggsave(tree_beta_plot, filename = here(plot_dir, "crumblr_cell_type_APOE_carrier_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE, low = "#398A84", high = "#D46B43") 
ggsave(forest_plot, filename = here(plot_dir, "crumblr_cell_type_APOE_carrier_forest.png"), width = 4, height = 7)


# combined_fig <- fig.vp +
#     theme(legend.position = "left") |
#     tree_beta_plot |
#     forest_plot
# 
# ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined_APOE_carrier.png"), width = 10)

## clr plot of top results
clr_boxplot_APOE_carrier_Excit  <- clr_prop_long |>
    filter(cell_type_anno %in% c("Excit.L2", "Excit.L6b")) |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~cell_type_anno) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier_Excit , filename = here(plot_dir, "clr_boxplot_APOE_carrier_Excit.png"), height = 4, width = 4)

clr_boxplot_APOE_carrier_OligoOPC <- clr_prop_long |>
    filter(cell_type_anno %in% c("Oligo.1", "Oligo.2", "Oligo.3", "OPC.5")) |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier_OligoOPC, filename = here(plot_dir, "clr_boxplot_APOE_carrier_OligoOPC.png"), height = 4, width = 8)

clr_boxplot_APOE_carrier_OligoOPC_color <- clr_prop_long |>
    filter(cell_type_anno %in% c("Oligo.1", "Oligo.2", "Oligo.3", "OPC.5")) |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1, aes(color = APOE)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    scale_color_manual(values = APOE_genotype_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier_OligoOPC, filename = here(plot_dir, "clr_boxplot_APOE_carrier_OligoOPC.png"), height = 4, width = 8)

prop_boxplot_APOE_carrier_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = APOE_carrier, y = prop, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1, aes(color = APOE)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    scale_color_manual(values = APOE_genotype_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(prop_boxplot_APOE_carrier_Oligo, filename = here(plot_dir, "prop_boxplot_APOE_carrier_Oligo.png"), height = 4, width = 8)

## Oligo3 vs. age

# clr_v_prop <- clr_prop_long |>
#     filter(cell_type_anno =="Oligo.3") |>
#     ggplot(aes(Age, CLR, color = APOE_carrier, shape = Sex)) +
#     geom_point()


#### fit exp_round ####
# Extract results for each cell type
topTable(fit, coef = "exp_roundround2", number = Inf)

#                         logFC     AveExpr            t     P.Value adj.P.Val         B
# Astro.5           1.106226545 -0.48566837  2.870510105 0.008739602 0.3321049 -2.623310
# Astro.1          -0.913041885  2.03579947 -1.891486214 0.071449543 0.7555904 -4.010021
# Astro.4          -1.171855941 -0.36993155 -1.873555972 0.073986158 0.7555904 -4.100112
# Oligo.5           0.695186652  1.16463328  1.678979887 0.106936170 0.7555904 -4.248310

topTable(fit, coef = "SexM", number = Inf)

# write.csv(diff_prop_APOE_carrier, file = here(data_dir, "diff_prop_APOE_carrier.csv"))


#### fit APOE ####
# 
# fit2 <- dream(cobj, ~ APOE + Age + Anc_Afr + exp_round , erc_info)
# fit2 <- eBayes(fit2)
# 
# # Test APOE carrier
# (geno <- colnames(fit2$design)[grep("APOE", colnames(fit2$design))])
# (diff_prop_APOE <- topTable(fit2, coef = geno, number = Inf))
# #                   APOEE2.E3    APOEE3.E4    APOEE4.E4    AveExpr         F    P.Value adj.P.Val
# # Excit.L6b        -0.3774346 -1.642397391 -1.017496429 -1.0645145 4.7032782 0.01140097 0.2374917
# # Oligo.2          -0.3884539 -1.142631475 -0.882662486  2.2535239 3.7707442 0.02588605 0.2374917
# # Excit.L2          0.1743530  0.984704658  0.653632584 -0.1314255 3.2310558 0.04272876 0.2374917
# 
# write.csv(diff_prop_APOE, file = here(data_dir, "diff_prop_APOE.csv"))
# 
# 
# # Perform multivariate test across the hierarchy
# topTable(fit, coef = "APOEE4/E4", number = Inf)
# 
# res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOEE4/E4")
# 
# # Plot hierarchy and testing results
# plotTreeTest(res)
# 
# # Plot hierarchy and regression coefficients
# plotTreeTestBeta(res)
# plotForest(res, hide = FALSE) 
# 
# combined_fig <- fig.vp +
#     theme(legend.position = "left") |
#     plotTreeTestBeta(res) +
#     theme(legend.position = "bottom", legend.box = "vertical") |
#     plotForest(res, hide = FALSE) 
# 
# ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined_APOEE4.E4.png"))
# 
# clr_boxplot_APOE_Excit.L6b <- clr_prop_long |>
#     filter(cell_type_anno == "Excit.L6b") |>
#     ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(aes(color = error), width = .1) +
#     facet_wrap(~cell_type_anno) +
#     scale_fill_manual(values = APOE_genotype_colors) +
#     theme_bw()
# 
# ggsave(clr_boxplot_APOE_Excit.L6b, filename = "clr_boxplot_APOE_Excit.L6b.png")
# 
# prop_boxplot_APOE_Excit.L6b <- clr_prop_long |>
#     filter(cell_type_anno == "Excit.L6b") |>
#     ggplot(aes(x = APOE, y = prop, fill = APOE)) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(aes(color = error), width = .1) +
#     facet_wrap(~cell_type_anno) +
#     scale_fill_manual(values = APOE_genotype_colors) +
#     theme_bw()
# 
# ggsave(prop_boxplot_APOE_Excit.L6b, filename = "prop_boxplot_APOE_Excit.L6b.png")

#### fit Sex ####

# Extract results for each cell type
(diff_prop_Sex <- topTable(fit, coef = "SexM", number = Inf))

#                        logFC     AveExpr           t     P.Value  adj.P.Val         B
# Oligo.5          -0.94666061  1.16152659 -3.29523642 0.003352853 0.07350752 -1.705243
# Oligo.3          -1.19927385  1.73901614 -3.22491808 0.003955976 0.07350752 -1.849009
# Micro.3           1.09017877  0.07707637  2.98973855 0.006836606 0.07350752 -2.419935
# Micro.4           0.85464395  0.09432697  2.93586040 0.007737634 0.07350752 -2.538259
# Astro.5          -0.89531742 -0.49334478 -2.76536537 0.011398815 0.08663099 -2.734235

write.csv(diff_prop_Sex, file = here(data_dir, "diff_prop_Sex.csv"))

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "SexM")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "sn_diff_prop_tree_test_Sex.csv"))

# Plot hierarchy and testing results
tree_fdr_plot <- plotTreeTest(res)
ggsave(tree_fdr_plot, filename = here(plot_dir, "crumblr_cell_type_Sex_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
tree_beta_plot <- plotTreeTestBeta(res, low = sex_colors[['F']], high = sex_colors[['M']]) +
    theme(legend.position = "left" 
          # legend.box = "vertical"
    )
ggsave(tree_beta_plot, filename = here(plot_dir, "crumblr_cell_type_Sex_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE, low = sex_colors[['F']], high = sex_colors[['M']]) 
ggsave(forest_plot, filename = here(plot_dir, "crumblr_cell_type_Sex_forest.png"), width = 4, height = 7)


# combined_fig <- fig.vp +
#     theme(legend.position = "left") |
#     tree_beta_plot |
#     forest_plot
# 
# ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined_Sex.png"), width = 10)

## clr plot of top results
clr_boxplot_Sex_Micro <- clr_prop_long |>
    filter(grepl("Micro", cell_type_anno)) |>
    ggplot(aes(x = Sex, y = CLR, fill = Sex)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = sex_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_Sex_Micro , filename = here(plot_dir, "clr_boxplot_Sex_Micro.png"), height = 4, width = 6)

clr_boxplot_Sex_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = Sex, y = CLR, fill = Sex)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_Sex_Oligo, filename = here(plot_dir, "clr_boxplot_Sex_Oligo.png"), height = 4, width = 6)

clr_boxplot_Sex_carrier_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    mutate(Sex_APOE = paste(Sex, APOE_carrier)) |>
    ggplot(aes(x = Sex_APOE, y = CLR, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Sex_carrier_Oligo, filename = here(plot_dir, "clr_boxplot_Sex_carrier_Oligo.png"), height = 4, width = 6)

clr_boxplot_Sex_carrier_Astro <- clr_prop_long |>
    filter(grepl("Astro", cell_type_anno)) |>
    mutate(Sex_APOE = paste(Sex, APOE_carrier)) |>
    ggplot(aes(x = Sex_APOE, y = CLR, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Sex_carrier_Astro, filename = here(plot_dir, "clr_boxplot_Sex_carrier_Astro.png"), height = 4, width = 6)

clr_boxplot_Sex_carrier_Micro <- clr_prop_long |>
    filter(grepl("Micro", cell_type_anno)) |>
    mutate(Sex_APOE = paste(Sex, APOE_carrier)) |>
    ggplot(aes(x = Sex_APOE, y = CLR, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Sex_carrier_Micro, filename = here(plot_dir, "clr_boxplot_Sex_carrier_Micro.png"), height = 4, width = 6)

prop_boxplot_Sex_carrier_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    mutate(Sex_APOE = paste(Sex, APOE_carrier)) |>
    ggplot(aes(x = Sex_APOE, y = prop, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(prop_boxplot_Sex_carrier_Oligo, filename = here(plot_dir, "prop_boxplot_Sex_carrier_Oligo.png"), height = 4, width = 6)

#### fit Age ####

# Extract results for each cell type
(diff_prop_Age<- topTable(fit, coef = "Age", number = Inf))

#                           logFC     AveExpr           t      P.Value  adj.P.Val          B
# Oligo.5           0.0445079005  1.16152659  3.88315362 0.0008221231 0.03124068 -0.6656332
# Oligo.3           0.0367116103  1.73901614  2.36749646 0.0272914275 0.48068363 -3.8912433
# Oligo.4           0.0317998131  1.65924744  2.21071638 0.0379487078 0.48068363 -4.1927394

write.csv(diff_prop_Age, file = here(data_dir, "diff_prop_Age.csv"))

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "Age")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "sn_diff_prop_tree_test_Age.csv"))

# Plot hierarchy and testing results
tree_fdr_plot <- plotTreeTest(res)
ggsave(tree_fdr_plot, filename = here(plot_dir, "crumblr_cell_type_Age_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
tree_beta_plot <- plotTreeTestBeta(res) +
    theme(legend.position = "left" 
          # legend.box = "vertical"
    )
ggsave(tree_beta_plot, filename = here(plot_dir, "crumblr_cell_type_Age_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE) 
ggsave(forest_plot, filename = here(plot_dir, "crumblr_cell_type_Age_forest.png"), width = 4, height = 7)

## clr plot of top results

clr_sactter_Age_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = Age, y = CLR)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~cell_type_anno, nrow = 1, scales = "free_y") +
    theme_bw() 

ggsave(clr_sactter_Age_Oligo, filename = here(plot_dir, "clr_sactter_Age_Oligo.png"), height = 4, width =10)

clr_sactter_Age_Oligo_APOE <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = Age, y = CLR, color = APOE_carrier)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~cell_type_anno, nrow = 1, scales = "free_y") +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() 

ggsave(clr_sactter_Age_Oligo_APOE, filename = here(plot_dir, "clr_sactter_Age_Oligo_APOE.png"), height = 4, width =10)

#### Run DREAM + eBayes on cobj interaction ####
## not enought data
# intertaction_fit <- dream(cobj, ~ APOE_carrier*Sex + Age + Anc_Afr + exp_round , erc_info) 
# intertaction_fit <- eBayes(fit = intertaction_fit)
# 
# (diff_prop_APOE_Sex <- topTable(intertaction_fit, coef = "APOE_carrierE4+:SexM", number = Inf))


# slurmjobs::job_single('22_crumblr_sn', create_shell = TRUE, memory = '25G', command = "Rscript 22_crumblr_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
