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

data_dir <- here("processed-data", "21_Xenium", "10_xenium_cell_types")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "10_xenium_cell_types")
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


table(rctd_data@results$results_df$spot_class)
# reject           singlet   doublet_certain doublet_uncertain 
# 5609            304401            122264             18033

table(rctd_data@results$results_df$first_type)

# Astro.1          Astro.2          Astro.3          Astro.4          Astro.5            Macro          Micro.1 
# 15125            18158              945             8954            42636             6285              291 
# Micro.2          Micro.3          Micro.4          Micro.5            OPC.1            OPC.2            OPC.3 
# 14780             2412              345            18178              120            11871               34 
# OPC.4            OPC.5          Oligo.1          Oligo.2          Oligo.3          Oligo.4          Oligo.5 
# 3814            11465            16888            11961            63536            37989             3599 
# Vasc.Endo          Vasc.PC        Vasc.VLMC         Excit.L2     Excit.L2_5.1     Excit.L2_5.2       Excit.L5.1 
# 36177            27006             8918            22124            10358             6601             8213 
# Excit.L5.2    Excit.L5_6_NP      Excit.L6_CT        Excit.L6b Inhib.Chandelier Inhib.Lamp5_Lhx6       Inhib.Pax6 
# 2341              983             5039             4974             1447             5284             1713 
# Inhib.Pvalb        Inhib.Sst        Inhib.Vip 
# 5247             4651             9845 

## spot class summary
rcdt_results_class_summary <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    count(sample_id, spot_class) |>
    group_by(sample_id) |>
    mutate(prop = n/sum(n))

## spot class summary plots
rctd_class_bar_plot <- rcdt_results_class_summary |>
    ggplot(aes(x = sample_id, y = n, fill = spot_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_bar_plot, filename = here(plot_dir, "xenium_rctd_class_bar_plot.png"))

rctd_class_prop_bar_plot <- rcdt_results_class_summary |>
    ggplot(aes(x = sample_id, y = prop, fill = spot_class)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(rctd_class_prop_bar_plot, filename = here(plot_dir, "xenium_rctd_class_prop_bar_plot.png"))

## cell type summary 
rcdt_results_singlets_summary <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    filter(spot_class == "singlet") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    group_by(sample_id, cell_type_anno = first_type) |>
    summarise(xenium_n = n()) |>
    group_by(sample_id) |>
    mutate(xenium_singlet_prop = xenium_n/sum(xenium_n))


rcdt_results_doublets <- rctd_data@results$results_df |> 
    rownames_to_column("barcode") |>
    filter(spot_class == "doublet_certain") |>
    mutate(sample_id = gsub("_.*", "", barcode)) |>
    select(barcode, sample_id, first_type, second_type) |>
    pivot_longer(!c(barcode, sample_id), names_to = "type", values_to = "cell_type_anno")

doublet_combos <- rcdt_results_doublets |> 
    group_by(barcode, sample_id) |> 
    summarise(doublet = paste0(sort(cell_type_anno), collapse = "-")) |>
    ungroup() |>
    count(doublet)

doublet_combos |> filter(grepl("Oligo.3", doublet)) |> arrange(-n)

# doublet                  n
# <chr>                <int>
# 1 Astro.5-Oligo.3       1714
# 2 Oligo.3-Excit.L6_CT   1615
# 3 Oligo.3-Excit.L6b     1608
# 4 Oligo.3-Excit.L2_5.1  1546
# 5 Oligo.3-Excit.L2_5.2  1022

rcdt_results_doublets_summary <- rcdt_results_doublets |>    
    group_by(sample_id, cell_type_anno) |>
    summarise(xenium_n = n())  |>
    group_by(sample_id) |>
    mutate(xenium_doublet_prop = xenium_n/sum(xenium_n))

rcdt_results_summary <- rcdt_results_singlets_summary |> 
    bind_rows(rcdt_results_doublets_summary)  |>
    group_by(sample_id, cell_type_anno) |>
    summarise(xenium_n = sum(xenium_n)) |>
    group_by(sample_id) |>
    mutate(xenium_prop = xenium_n/sum(xenium_n))


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

#### PCA ####
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

mod <- model.matrix( ~ APOE_carrier + Sex + Age + Anc_Afr + exp_round , erc_info)
cleanY_cobj <- cleaningY(cobj$E , mod, P=2)
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


