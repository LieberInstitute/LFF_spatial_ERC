## Louise Huuki-Myers, May 2026
## Explore RCTD results

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")

library("tidyverse")
library("crumblr")
library("variancePartition")

data_dir <- here("processed-data", "21_Xenium", "11_xenium_crumblr")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "11_xenium_crumblr")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))

message(Sys.time(), "- Load rctd data")
rctd_data <- qs_read(here("processed-data", "21_Xenium", "09_xenium_label_transfer_RCTD","rctd_results_xenium.qs2"))

names(rctd_data@results)

head(rctd_data@results$results_df)

# spot_class first_type second_type first_class second_class min_score singlet_score conv_all
# Br1039_aaaaabke-1 doublet_certain  Vasc.Endo     Oligo.4       FALSE         TRUE 181.70462      215.9917     TRUE
# Br1039_aaaaepee-1         singlet    Astro.4     Oligo.3        TRUE        FALSE 117.32119      133.8569     TRUE
# Br1039_aaagdfnc-1         singlet    Astro.2     Oligo.2       FALSE        FALSE 235.05667      237.8763     TRUE
# Br1039_aaajpooh-1         singlet    Astro.5     Oligo.3       FALSE        FALSE  94.88785      107.9397     TRUE
# Br1039_aaalplmo-1          reject    Astro.5     Oligo.3       FALSE        FALSE 157.67283      197.6728     TRUE
# Br1039_aaamomng-1 doublet_certain    Oligo.3     Micro.2       FALSE        FALSE 194.25099      226.9972     TRUE
# conv_doublet
# Br1039_aaaaabke-1         TRUE
# Br1039_aaaaepee-1         TRUE
# Br1039_aaagdfnc-1         TRUE
# Br1039_aaajpooh-1         TRUE
# Br1039_aaalplmo-1         TRUE
# Br1039_aaamomng-1         TRUE

# boolean variables `first_class` and
# `second_class` represent, for the first and second cell type,
# respectively, whether an assignment was made on the class level



#### Compare proportions to sn-RNA seq data ####

# cell_type_proportions
load(here("processed-data", "04_snRNA-seq", "28_subcluster_update_sce","cell_type_proportions.Rdata"), verbose = TRUE)

cell_type_proportions_xenium <- cell_type_proportions |>
    mutate(cell_type_broad = factor(gsub("\\..*", "", cell_type_anno), levels = names(cell_type_colors$broad)), .before = 10) |>
    select(-exp_round, -seq_round) |> 
    inner_join(rcdt_results_summary)

sn_vs_xenium_n_scatter <- cell_type_proportions_xenium |>
    ggplot(aes(n, xenium_n, color = cell_type_anno)) +
    geom_point() +
    scale_color_manual(values = cell_type_colors$anno) + 
    facet_wrap(~cell_type_broad, scales = "free") +
    theme_bw()

ggsave(sn_vs_xenium_n_scatter, filename = here(plot_dir, "xenium_rctd_vs_sn_cell_type_n.png"), height = 10, width = 10)


sn_vs_xenium_prop_scatter <- cell_type_proportions_xenium |>
    ggplot(aes(prop, xenium_prop, color = cell_type_anno)) +
    geom_point() +
    scale_color_manual(values = cell_type_colors$anno) + 
    facet_wrap(~cell_type_broad, scales = "free") +
    theme_bw()

ggsave(sn_vs_xenium_prop_scatter, filename = here(plot_dir, "xenium_rctd_vs_sn_cell_type_prop.png"), height = 10, width = 10)

sn_vs_xenium_prop_boxplot <- cell_type_proportions_xenium |>
    ungroup() |>
    select(sample_id, cell_type_anno, cell_type_broad, sn_prop = prop, xenium_prop) |>
    pivot_longer(!c(sample_id, cell_type_anno, cell_type_broad), names_to = "assay", values_to = "prop") |>
    ggplot(aes(x = cell_type_anno, y = prop, color = assay)) +
    geom_boxplot() +
    facet_wrap(~cell_type_broad, scales = "free", space = "free_x")  +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sn_vs_xenium_prop_boxplot, filename = here(plot_dir, "xenium_rctd_vs_sn_cell_type_prop_boxplot.png"), height = 7, width = 15)

#### CRUMBLR cell type proportions ####

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- cell_type_proportions_xenium |> #already grouped by sample
    mutate(error = sqrt((xenium_prop * (1-xenium_prop))/sum(xenium_n)),
           anno = paste0(xenium_n,"/", sum(xenium_n))) |>
    select(sample_id, cell_type_anno, prop = xenium_prop, error, anno)

cell_counts <- cell_type_proportions_xenium |> 
    ungroup() |>
    select(sample_id, cell_type_anno, xenium_n) |>
    pivot_wider(names_from = cell_type_anno, values_from = xenium_n)  |>
    column_to_rownames("sample_id")

cell_counts[is.na(cell_counts)] <- 0

## create info matrix
xenium_experiment_info <- read.csv(here("processed-data", "21_Xenium", "07_xenium_build_spe", "xenium_experiment_details.csv"))

erc_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv")) |>
    select(sample_id, Age, Sex, Ancestry, Anc_Afr, APOE_carrier, APOE) |> 
    inner_join(xenium_experiment_info |> select(sample_id = BrNum, Run, chip)) |>
    filter(sample_id %in% rownames(cell_counts)) |>
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
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | Run) + (1 | chip)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)
fig.vp

ggsave(fig.vp + 
           theme_bw() + 
           theme(legend.position = "bottom"),
       filename = here(plot_dir, "xenium_crumblr_cell_type_vp.png"), width = 4.5, height = 7)

#### variable correlation ####
C <- canCorPairs(form, erc_info)
# APOE and APOE_carrier
pdf(here(plot_dir, "xenium_variable_CorrMatrix.pdf"), height = 5, width = 5)
plotCorrMatrix(C)
dev.off()

#### Crumblr PCA ####
pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
#   shape by Stimulated vs unstimulated
cobj_pca <- ggplot(df_pca, aes(PC1, PC2, color = Run, shape = APOE)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

ggsave(cobj_pca, filename = here(plot_dir, "xenium_crumblr_pca.png"))

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
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(clr_boxplot_APOE, filename = here(plot_dir, "xenium_clr_boxplot_APOE.png"), width = 10)

clr_boxplot_APOE_jitter <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_anno) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(clr_boxplot_APOE_jitter, filename = here(plot_dir, "xenium_clr_boxplot_APOE_jitter.png"), width = 10)

clr_boxplot_APOE_carrier <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_anno) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_carrier, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier.png"))

clr_boxplot_chip <- clr_prop_long |>
    ggplot(aes(x = chip, y = CLR, fill = Run)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cell_type_anno) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_chip, filename = here(plot_dir, "xenium_clr_boxplot_chip.png"), width = 10)

#### Run DREAM + eBayes on cobj ####
fit <- dream(cobj, ~ APOE_carrier + Sex + Age + Anc_Afr + chip , erc_info) 
fit <- eBayes(fit = fit)

## Load Hierarchical clustering
load(here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster","sn_subcluster_hierarchical_cluster.Rdata"), verbose = TRUE)

## fit APOE_carrier 
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance

# Extract results for each cell type
(diff_prop_APOE_carrier <- topTable(fit, coef = "APOE_carrierE4+", number = Inf))

#                        logFC      AveExpr           t     P.Value  adj.P.Val         B
# Micro.3          -1.13696090 -0.565597771 -3.78059401 0.003638408 0.06372119 -1.577516
# Micro.4          -0.89867263 -2.188224745 -3.77484929 0.003672728 0.06372119 -1.654127
# Micro.1          -0.73505706 -2.108237036 -3.58366692 0.005030621 0.06372119 -1.862773
# Excit.L2          0.75563619  1.523895447  2.66028590 0.024006653 0.22806320 -3.201613

write.csv(diff_prop_APOE_carrier, file = here(data_dir, "xenium_diff_prop_APOE_carrier.csv"))

#### cleanY on cobj  APOE_carrier ####

mod <- model.matrix( ~ APOE_carrier + Sex + Age + Anc_Afr + chip , erc_info)
cleanY_cobj <- jaffelab::cleaningY(cobj$E , mod, P=2)
head(cleanY_cobj)

clean_clr_prop_long <- cleanY_cobj |>
    as.data.frame() |>
    rownames_to_column("cell_type_anno") |>
    pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR_cleanY")

clr_prop_long <- clr_prop_long |> left_join(clean_clr_prop_long)

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "xenium_diff_prop_tree_test_APOE_carrier.csv"))

# Plot hierarchy and testing results
tree_fdr_plot <- plotTreeTest(res)
ggsave(tree_fdr_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_APOE_carrier_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
# TODO fdr.cutoff = 0.1
tree_beta_plot <- plotTreeTestBeta(res, low = "#398A84", high = "#D46B43") +
    theme(legend.position = "left" 
          # legend.box = "vertical"
    )

ggsave(tree_beta_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_APOE_carrier_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE, low = "#398A84", high = "#D46B43") 
ggsave(forest_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_APOE_carrier_forest.png"), width = 4, height = 7)

## clr plot of top results
clr_boxplot_APOE_carrier_Excit  <- clr_prop_long |>
    filter(cell_type_anno %in% c("Excit.L2", "Excit.L2_5.1")) |>
    ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~cell_type_anno, scales = "free_y") +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier_Excit , filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_Excit.png"), height = 4, width = 4)

clr_boxplot_APOE_carrier_OligoOPC <- clr_prop_long |>
    filter(cell_type_anno %in% c("Oligo.1", "Oligo.2", "Oligo.3", "OPC.5")) |>
    ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_APOE_carrier_OligoOPC, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_OligoOPC.png"), height = 4, width = 8)

clr_boxplot_APOE_carrier_OligoOPC_color <- clr_prop_long |>
    filter(cell_type_anno %in% c("Oligo.1", "Oligo.2", "Oligo.3", "OPC.5")) |>
    ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
    # geom_boxplot() +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    geom_jitter(width = .1, aes(color = APOE)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = APOE_carrier_colors) +
    scale_color_manual(values = APOE_genotype_colors) +
    theme_bw() 

ggsave(clr_boxplot_APOE_carrier_OligoOPC_color, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_OligoOPC_color.png"), height = 4, width = 8)

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

ggsave(prop_boxplot_APOE_carrier_Oligo, filename = here(plot_dir, "xenium_prop_boxplot_APOE_carrier_Oligo.png"), height = 4, width = 8)

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


#### Plot cleanY on cobj Sex ####

mod <- model.matrix( ~ APOE_carrier + Sex +  Age + Anc_Afr + exp_round , erc_info)
cleanY_cobj_sex <- cleaningY(cobj$E , mod, P=3)

clr_prop_long <- clr_prop_long |> 
    left_join(cleanY_cobj_sex |>
                  as.data.frame() |>
                  rownames_to_column("cell_type_anno") |>
                  pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR_cleanY_Sex")) 

## clr plot of top results
clr_boxplot_Sex_Micro <- clr_prop_long |>
    filter(grepl("Micro", cell_type_anno)) |>
    ggplot(aes(x = Sex, y = CLR_cleanY_Sex, fill = Sex)) +
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
    ggplot(aes(x = Sex, y = CLR_cleanY_Sex, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = sex_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(clr_boxplot_Sex_Oligo, filename = here(plot_dir, "clr_boxplot_Sex_Oligo.png"), height = 4, width = 6)


#### Plot cleanY on cobj Sex + APOE carrier ####

mod <- model.matrix( ~ Sex + APOE_carrier + Age + Anc_Afr + exp_round , erc_info)
cleanY_cobj_carrier_sex <- cleaningY(cobj$E , mod, P=3)

clr_prop_long <- clr_prop_long |> 
    left_join(cleanY_cobj_sex |>
                  as.data.frame() |>
                  rownames_to_column("cell_type_anno") |>
                  pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR_cleanY_carrier_Sex")) |>
    mutate(Sex_APOE = paste(Sex, APOE_carrier))

## plots 
clr_boxplot_Sex_carrier_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = Sex_APOE, y = CLR_cleanY_carrier_Sex, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Sex_carrier_Oligo, filename = here(plot_dir, "clr_boxplot_Sex_carrier_Oligo.png"), height = 4, width = 6)

## plot prop
prop_boxplot_Sex_carrier_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
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


clr_boxplot_Sex_carrier_Astro <- clr_prop_long |>
    filter(grepl("Astro", cell_type_anno)) |>
    ggplot(aes(x = Sex_APOE, y = CLR_cleanY_carrier_Sex, fill = Sex)) +
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
    ggplot(aes(x = Sex_APOE, y = CLR_cleanY_carrier_Sex, fill = Sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .1, aes(color = APOE_carrier)) +
    facet_wrap(~cell_type_anno, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_Sex_carrier_Micro, filename = here(plot_dir, "clr_boxplot_Sex_carrier_Micro.png"), height = 4, width = 6)


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

#### Plot cleanY CLR for Age ####

mod <- model.matrix( ~ Age + APOE_carrier + Sex + Anc_Afr + exp_round , erc_info)
cleanY_cobj_Age <- cleaningY(cobj$E , mod, P=2)

clr_prop_long$CLR_cleanY_Age <- NULL

clr_prop_long <- clr_prop_long |> 
    left_join(cleanY_cobj_Age |>
                  as.data.frame() |>
                  rownames_to_column("cell_type_anno") |>
                  pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR_cleanY_Age")) 


clr_sactter_Age_Oligo <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = Age, y = CLR_cleanY_Age)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~cell_type_anno, nrow = 1) +
    theme_bw() 

ggsave(clr_sactter_Age_Oligo, filename = here(plot_dir, "clr_sactter_Age_Oligo.png"), height = 4, width =10)


#### Plot cleanY CLR for Age & APOE ####

# mod <- model.matrix( ~ Age + APOE_carrier + Sex + Anc_Afr + exp_round , erc_info)
cleanY_cobj_Age <- cleaningY(cobj$E , mod, P=3)

clr_prop_long$CLR_cleanY_carrier_Age <- NULL

clr_prop_long <- clr_prop_long |> 
    left_join(cleanY_cobj_Age |>
                  as.data.frame() |>
                  rownames_to_column("cell_type_anno") |>
                  pivot_longer(!cell_type_anno, names_to = "sample_id", values_to = "CLR_cleanY_carrier_Age")) 


clr_sactter_Age_Oligo_APOE <- clr_prop_long |>
    filter(grepl("Oligo", cell_type_anno)) |>
    ggplot(aes(x = Age, y = CLR_cleanY_carrier_Age, color = APOE_carrier)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~cell_type_anno, nrow = 1) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave(clr_sactter_Age_Oligo_APOE, filename = here(plot_dir, "clr_sactter_Age_Oligo_APOE.png"), height = 4, width =10)


# slurmjobs::job_single('22_crumblr_sn', create_shell = TRUE, memory = '25G', command = "Rscript 22_crumblr_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

