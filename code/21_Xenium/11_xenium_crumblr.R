## Louise Huuki-Myers, May 2026
## Explore RCTD results
## Crumblr proportion tests on ERC Xenium data, by cell type or xSpD
## Usage: Rscript 11_xenium_crumblr.R --cluster cell_type   (default)
##        Rscript 11_xenium_crumblr.R --cluster xSpD

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("jaffelab")
library("tidyverse")
library("crumblr")
library("variancePartition")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("cluster", "c", "1", "character", "Name of cluster: cell_type (default) or xSpD"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
if(is.null(opt$cluster)) opt$cluster <- "cell_type"
stopifnot(opt$cluster %in% c("cell_type", "xSpD"))

is_cell_type <- opt$cluster == "cell_type"
clus_var <- ifelse(is_cell_type, "cell_type_anno", "xSpD")
clus_suffix <- ifelse(is_cell_type, "", paste0("-", opt$cluster))

## Hierarchical clustering used as the tree for treeTest
tree_file <- switch(opt$cluster,
    cell_type = here("processed-data", "04_snRNA-seq", "32_sn_subcluster_hierarchical_cluster", "sn_subcluster_hierarchical_cluster.Rdata"),
    xSpD = here("processed-data", "21_Xenium", "24_xenium_xSpD_pseudobulk_hclust", "xSpD_hierarchical_cluster.Rdata")
)
stopifnot(file.exists(tree_file))

data_dir <- here("processed-data", "21_Xenium", paste0("11_xenium_crumblr", clus_suffix))
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", paste0("11_xenium_crumblr", clus_suffix))
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### Load data ####
message(Sys.time(), "- Load xenium data")
## xSpD spe is not filtered to singlets (see 14_xenium_pseudobulk_model_register.R)
spe_fn <- switch(opt$cluster,
    cell_type = here("processed-data", "21_Xenium", "10_xenium_cell_types", "spe_xenium_cell_types.qs2"),
    xSpD = here("processed-data", "21_Xenium", "13_xenium_bansky_embedding", "spe_xenium_bansky.qs2")
)
spe <- qs_read(spe_fn)

#### Xenium proportions per sample ####
message(Sys.time(), "- Summarize xenium proportions")
rcdt_results_summary <- as.data.frame(colData(spe)) |>
    group_by(sample_id, cluster = .data[[clus_var]]) |>
    summarise(xenium_n = n(), .groups = "drop_last") |>
    group_by(sample_id) |>
    mutate(xenium_prop = xenium_n/sum(xenium_n))

prop_xenium <- rcdt_results_summary

if(is_cell_type) {
    #### Compare proportions to sn-RNA seq data ####
    message(Sys.time(), "- Compare with snRNA prop")

    # cell_type_proportions
    load(here("processed-data", "04_snRNA-seq", "28_subcluster_update_sce","cell_type_proportions.Rdata"), verbose = TRUE)

    cell_type_proportions_xenium <- cell_type_proportions |>
        mutate(cell_type_broad = factor(gsub("\\..*", "", cell_type_anno), levels = names(cell_type_colors$broad)), .before = 10) |>
        select(-exp_round, -seq_round) |>
        inner_join(rcdt_results_summary |> rename(cell_type_anno = cluster))

    write_csv(cell_type_proportions_xenium, here(data_dir, "cell_type_proportions_xenium.csv"))

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

    ## only keep xenium cell types that have a snRNA match
    prop_xenium <- cell_type_proportions_xenium |>
        rename(cluster = cell_type_anno)
}

#### CRUMBLR proportions ####

## calc error formula for standard error of propotion SE = sqrt(p*(1-p)/n)
error <- prop_xenium |>
    group_by(sample_id) |>
    mutate(error = sqrt((xenium_prop * (1-xenium_prop))/sum(xenium_n)),
           anno = paste0(xenium_n,"/", sum(xenium_n))) |>
    ungroup() |>
    select(sample_id, cluster, prop = xenium_prop, error, anno)

cell_counts <- prop_xenium |> 
    ungroup() |>
    select(sample_id, cluster, xenium_n) |>
    pivot_wider(names_from = cluster, values_from = xenium_n)  |>
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
    rownames_to_column("cluster") |>
    pivot_longer(!cluster, names_to = "sample_id", values_to = "CLR") |>
    left_join(erc_info |> rownames_to_column("sample_id")) |>
    left_join(error) |>
    replace_na(list(prop = 0, error = 0))

# Partition variance into components for sample
form <- ~ (1 | APOE_carrier) + (1 | Sex) + Age + Anc_Afr + (1 | Run) + (1 | chip)
vp <- fitExtractVarPartModel(cobj, form, erc_info)

# Plot variance fractions
fig.vp <- plotPercentBars(vp)

ggsave(fig.vp + 
           theme_bw() + 
           theme(legend.position = "bottom"),
       filename = here(plot_dir, "xenium_crumblr_cell_type_vp.png"), width = 4.5, height = 7)

#### variable correlation ####
message(Sys.time(), "- Check variable cor")

C <- canCorPairs(form, erc_info)
# APOE and APOE_carrier
pdf(here(plot_dir, "xenium_variable_CorrMatrix.pdf"), height = 5, width = 5)
plotCorrMatrix(C)
dev.off()


# High colinearity between variables:
#     Run and chip

#### Crumblr PCA ####
message(Sys.time(), "- CRUMBLR PCA")

pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, erc_info, by = "row.names")

# Plot PCA
cobj_pca <- ggplot(df_pca, aes(PC1, PC2, color = Run, shape = APOE)) +
    geom_point(size = 3) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    scale_color_discrete(name = "Subject") +
    xlab("PC1") +
    ylab("PC2")

ggsave(cobj_pca, filename = here(plot_dir, "xenium_crumblr_pca.png"))

#### CLR plots ####
message(Sys.time(), "- CLR plots")

if(is_cell_type) {
    clr_prop_long |> filter(cluster == "Oligo.3") |> arrange(-prop)
}

# clr_v_prop <- clr_prop_long |>
#     ggplot(aes(prop, CLR, color = error)) +
#     geom_point()

clr_boxplot_APOE <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot() +
    facet_wrap(~cluster) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(clr_boxplot_APOE, filename = here(plot_dir, "xenium_clr_boxplot_APOE.png"), width = 10)

clr_boxplot_APOE_jitter <- clr_prop_long |>
    ggplot(aes(x = APOE, y = CLR, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cluster) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(clr_boxplot_APOE_jitter, filename = here(plot_dir, "xenium_clr_boxplot_APOE_jitter.png"), width = 10)

clr_boxplot_APOE_carrier <- clr_prop_long |>
    ggplot(aes(x = APOE_carrier, y = CLR, fill = APOE_carrier)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cluster) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()

ggsave(clr_boxplot_APOE_carrier, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier.png"))

clr_boxplot_chip <- clr_prop_long |>
    ggplot(aes(x = chip, y = CLR, fill = Run)) +
    geom_boxplot() +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter(aes(color = error), width = .1) +
    facet_wrap(~cluster) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(clr_boxplot_chip, filename = here(plot_dir, "xenium_clr_boxplot_chip.png"))

#### Save CLR + proprotion data ####
message(Sys.time(), "Save data")

write_csv(clr_prop_long, here(data_dir, "xenium_clr_prop_long.csv"))

#### Run DREAM + eBayes on cobj ####
message(Sys.time(), "- Run DREAM + eBayes")

fit <- dream(cobj, ~ APOE_carrier + Age + Anc_Afr + chip , erc_info) 
fit <- eBayes(fit = fit)

## Load hierarchical clustering
load(tree_file, verbose = TRUE)

## fit APOE_carrier 
# Perform regression on each cell type separately
#  then use eBayes to shrink residual variance

# Extract results for each cluster
(diff_prop_APOE_carrier <- topTable(fit, coef = "APOE_carrierE4+", number = Inf))

#                         logFC    AveExpr           t    P.Value adj.P.Val         B
# Astro.5           0.536476892  2.6131444  3.01008684 0.01196180 0.3026889 -3.403692
# Micro.3          -1.076747228 -1.0925486 -2.69203615 0.02108404 0.3026889 -3.302067
# OPC.5             0.451219869  1.3098030  2.44133565 0.03290517 0.3026889 -3.753504
# Micro.1          -1.269563402 -3.4182770 -2.37828394 0.03677932 0.3026889 -3.606209
# Excit.L2          0.679411326  1.5634345  2.17298967 0.05269257 0.3026889 -4.005760

write.csv(diff_prop_APOE_carrier, file = here(data_dir, "xenium_diff_prop_APOE_carrier.csv"))


#### cleanY on cobj  APOE_carrier ####
message(Sys.time(), "- Run cleanY")

mod <- model.matrix( ~ APOE_carrier + Age + Anc_Afr + chip , erc_info)
cleanY_cobj <- jaffelab::cleaningY(cobj$E , mod, P=2)
head(cleanY_cobj)

clean_clr_prop_long <- cleanY_cobj |>
    as.data.frame() |>
    rownames_to_column("cluster") |>
    pivot_longer(!cluster, names_to = "sample_id", values_to = "CLR_cleanY")

clr_prop_long <- clr_prop_long |> left_join(clean_clr_prop_long)

#### Fit APOE carrier ####
message('Fit APOE carrier')

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "APOE_carrierE4+")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "xenium_diff_prop_tree_test_APOE_carrier.csv"))

# Plot hierarchy and testing results
nominal_fdr_cut <- res_tb |> filter(pvalue < 0.1) |> slice_max(FDR) |> pull(FDR) |> max()
# nominal_fdr_cut <- nominal_fdr_cut + 0.0001 ## add buffer for plot cutoff

res_tb |> filter(FDR < nominal_fdr_cut) |> select(label, beta,pvalue, FDR)

res_tb |> count(FDR < nominal_fdr_cut)

tree_fdr_plot <- plotTreeTest(res, fdr.cutoff = nominal_fdr_cut)
ggsave(tree_fdr_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_APOE_carrier_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
tree_beta_plot <- plotTreeTestBeta(res, low = "#398A84", high = "#D46B43", fdr.cutoff = nominal_fdr_cut) +
    theme(legend.position = "left" 
          # legend.box = "vertical"
    )

ggsave(tree_beta_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_APOE_carrier_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE, low = "#398A84", high = "#D46B43") 
ggsave(forest_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_APOE_carrier_forest.png"), width = 4, height = 7)

if(is_cell_type) {
    ## clr plot of top results
    clr_boxplot_APOE_carrier_Excit  <- clr_prop_long |>
        filter(cluster %in% c("Excit.L2", "Excit.L2_5.1")) |>
        ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
        # geom_boxplot() +
        geom_boxplot(outlier.shape = NA) +
        # geom_jitter(aes(color = error), width = .1) +
        geom_jitter(width = .1) +
        facet_wrap(~cluster, scales = "free_y") +
        scale_fill_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "None") +
        labs(y = "Xenium CLR cleanY")

    ggsave(clr_boxplot_APOE_carrier_Excit , filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_Excit.png"), height = 4, width = 4)

    clr_boxplot_APOE_carrier_OligoOPC <- clr_prop_long |>
        filter(cluster %in% c("Oligo.1", "Oligo.2", "Oligo.3", "OPC.5")) |>
        ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
        # geom_boxplot() +
        geom_boxplot(outlier.shape = NA) +
        # geom_jitter(aes(color = error), width = .1) +
        geom_jitter(width = .1) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "None")  +
        labs(y = "Xenium CLR cleanY")

    ggsave(clr_boxplot_APOE_carrier_OligoOPC, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_OligoOPC.png"), height = 4, width = 8)

    clr_boxplot_APOE_carrier_OligoOPC_color <- clr_prop_long |>
        filter(cluster %in% c("Oligo.1", "Oligo.2", "Oligo.3", "OPC.5")) |>
        ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
        # geom_boxplot() +
        geom_boxplot(outlier.shape = NA) +
        # geom_jitter(aes(color = error), width = .1) +
        geom_jitter(width = .1, aes(color = APOE)) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = APOE_carrier_colors) +
        scale_color_manual(values = APOE_genotype_colors) +
        theme_bw() 

    ggsave(clr_boxplot_APOE_carrier_OligoOPC_color, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_OligoOPC_color.png"), height = 4, width = 8)

    prop_boxplot_APOE_carrier_Oligo <- clr_prop_long |>
        filter(grepl("Oligo", cluster)) |>
        ggplot(aes(x = APOE_carrier, y = prop, fill = APOE_carrier)) +
        # geom_boxplot() +
        geom_boxplot(outlier.shape = NA) +
        # geom_jitter(aes(color = error), width = .1) +
        geom_jitter(width = .1, aes(color = APOE)) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = APOE_carrier_colors) +
        scale_color_manual(values = APOE_genotype_colors) +
        theme_bw() +
        # theme(legend.position = "None") +
        labs(y = "Xenium prop")

    ggsave(prop_boxplot_APOE_carrier_Oligo, filename = here(plot_dir, "xenium_prop_boxplot_APOE_carrier_Oligo.png"), height = 4, width = 8)

    prop_boxplot_APOE_carrier_OPC <- clr_prop_long |>
        filter(grepl("OPC", cluster)) |>
        ggplot(aes(x = APOE_carrier, y = prop, fill = APOE_carrier)) +
        # geom_boxplot() +
        geom_boxplot(outlier.shape = NA) +
        # geom_jitter(aes(color = error), width = .1) +
        geom_jitter(width = .1, aes(color = APOE)) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = APOE_carrier_colors) +
        scale_color_manual(values = APOE_genotype_colors) +
        theme_bw() +
        # theme(legend.position = "None") +
        labs(y = "Xenium prop")

    ggsave(prop_boxplot_APOE_carrier_OPC, filename = here(plot_dir, "xenium_prop_boxplot_APOE_carrier_OPC.png"), height = 4, width = 8)

    clr_boxplot_APOE_carrier_Micro <- clr_prop_long |>
        filter(grepl("Micro", cluster)) |>
        ggplot(aes(x = APOE_carrier, y = CLR_cleanY, fill = APOE_carrier)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1, aes(color = APOE)) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = APOE_carrier_colors) +
        scale_color_manual(values = APOE_genotype_colors) +
        theme_bw() +
        theme(legend.position = "None")

    ggsave(clr_boxplot_APOE_carrier_Micro, filename = here(plot_dir, "xenium_clr_boxplot_APOE_carrier_Micro.png"), height = 4, width = 8)
}

if(is_cell_type) {
    #### compare to snRNA-seq crumblr Carrier tests ####
    message(Sys.time(), "- compare to snRNA-seq crumblr Carrier tests")

    sn_diff_prop_APOE_carrier <- read.csv(here("processed-data", "04_snRNA-seq", "22_crumblr_sn", "sn_diff_prop_tree_test_APOE_carrier.csv")) 

    x_name <- function(name) paste0("x_", name)
    p_cutoff <- 0.1

    diff_prop_tree_APOE_carrier_compare <- res_tb |> 
        rename_at(5:12, x_name) |>
        select(-branch.length) |> ## stored branch lengths cause buggy join
        full_join(sn_diff_prop_APOE_carrier) |>
        mutate(
            test = "APOE_carrier",
            label_multi = label |> str_split("/"),
            n_labels = lengths(label_multi),
            label_simple = map_chr(label_multi, ~ .x |>
                                       str_remove_all("\\.[A-Za-z0-9_]+") |>
                                       unique() |>
                                       paste(collapse = "/")
            ),
            label_simple2 = ifelse(n_labels < 3, label, paste0(label_simple, '[', node,']')),
            cell_type_broad = ifelse(grepl("/", label_simple), "MULTI", label_simple),
            signif_cat = case_when(
                pvalue < p_cutoff  & x_pvalue < p_cutoff ~ "both",
                pvalue < p_cutoff  & x_pvalue >= p_cutoff ~ "sn only",
                pvalue >= p_cutoff & x_pvalue < p_cutoff ~ "xenium only",
                pvalue >= p_cutoff & x_pvalue >= p_cutoff ~ "neither"),
            type = ifelse(n_labels > 1, "node", "leaf")
        )

    # diff_prop_tree_APOE_carrier_compare |> filter(n_labels == 2) |> select(label, n_labels, label_multi, label_simple, label_simple2, cell_type_broad) 
    # diff_prop_tree_APOE_carrier_compare |> dplyr::count(label_simple)
    # 
    # diff_prop_tree_APOE_carrier_compare$label

    signif_colors <- c(both = "purple", `xenium only` = "blue", `sn only` = "red")

    tree_compare_beta_scatter <- diff_prop_tree_APOE_carrier_compare |>
        ggplot(aes(beta, x_beta, color = signif_cat)) +
        geom_point() +
        geom_abline(color = "red", linetype = "dashed") +
        ggrepel::geom_text_repel(aes(label = label_simple2), size = 1.7) +
        scale_color_manual(values = signif_colors) +
        theme_bw() +
        labs(x = "snRNA-seq beta", y = "Xenium beta") 

    ggsave(tree_compare_beta_scatter, 
           filename = here(plot_dir, "tree_compare_beta_scatter_APOE_carrier.png"), width = 6, height = 5)

    ggsave(tree_compare_beta_scatter + facet_wrap(~type, nrow = 1), 
           filename = here(plot_dir, "tree_compare_beta_scatter_APOE_carrier_facet.png"), width = 11, height = 5)


    tree_compare_beta_scatter_ct_broad <- diff_prop_tree_APOE_carrier_compare |>
        ggplot(aes(beta, x_beta, color = cell_type_broad, shape = signif_cat)) +
        geom_point() +
        geom_abline(color = "red", linetype = "dashed") +
        ggrepel::geom_text_repel(aes(label = label_simple2), size = 1.7) +
        scale_color_manual(values = cell_type_colors$broad) +
        theme_bw()  +
        labs(x = "snRNA-seq beta", y = "Xenium beta") 

    ggsave(tree_compare_beta_scatter_ct_broad, filename = here(plot_dir, "tree_compare_beta_scatter_APOE_carrier_ct_broad.png"), width = 6, height = 5)


    tree_compare_stat_scatter_ct_broad <- diff_prop_tree_APOE_carrier_compare |>
        # filter(cell_type_broad != "MULTI") |>
        ggplot(aes(stat, x_stat, color = signif_cat, shape = pvalue < 0.1)) +
        geom_point() +
        # geom_abline(color = "red", linetype = "dashed") +
        ggrepel::geom_text_repel(aes(label = label_simple2), size = 2) +
        # scale_color_manual(values = cell_type_colors$broad) +
        theme_bw() +
        labs(x = "snRNA-seq stat", y = "Xenium stat") 
        # facet_wrap(~type)

    ggsave(tree_compare_stat_scatter_ct_broad, filename = here(plot_dir, "tree_compare_stat_scatter_APOE_carrier_ct_broad.png"))
}

#### fit Sex ####

# Extract results for each cluster
(diff_prop_Sex <- topTable(fit, coef = "SexM", number = Inf))

#                   logFC     AveExpr           t     P.Value  adj.P.Val         B
# Macro             1.99051694 -0.6335688  2.43415100 0.03533249 0.7314477 -4.509818
# OPC.1             1.42682280 -4.3314793  2.03835503 0.06900192 0.7314477 -4.551939
# Micro.2           0.84770877  1.3614673  1.99108878 0.07465394 0.7314477 -4.560956
# Micro.1           1.71162443 -3.4181275  1.81243561 0.10018125 0.7314477 -4.529114
# Excit.L2_5.1     -1.01489921  0.7277704 -1.77773943 0.10599298 0.7314477 -4.557097

write.csv(diff_prop_Sex, file = here(data_dir, "diff_prop_Sex.csv"))

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "SexM")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "xenium_diff_prop_tree_test_Sex.csv"))

# res_tb <- read_csv(here(data_dir, "xenium_diff_prop_tree_test_Sex.csv"))

# Plot hierarchy and testing results
tree_fdr_plot <- plotTreeTest(res)
ggsave(tree_fdr_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_Sex_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
tree_beta_plot <- plotTreeTestBeta(res, low = sex_colors[['F']], high = sex_colors[['M']]) +
    theme(legend.position = "left"
          # legend.box = "vertical"
    )
ggsave(tree_beta_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_Sex_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE, low = sex_colors[['F']], high = sex_colors[['M']])
ggsave(forest_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_Sex_forest.png"), width = 4, height = 7)


# combined_fig <- fig.vp +
#     theme(legend.position = "left") |
#     tree_beta_plot |
#     forest_plot
#
# ggsave(combined_fig, filename = here(plot_dir, "crumblr_cell_type_combined_Sex.png"), width = 10)


#### Plot cleanY on cobj Sex ####

mod <- model.matrix( ~ APOE_carrier + Sex +  Age + Anc_Afr + chip , erc_info)
cleanY_cobj_sex <- jaffelab::cleaningY(cobj$E , mod, P=3)

clr_prop_long <- clr_prop_long |>
    left_join(cleanY_cobj_sex |>
                  as.data.frame() |>
                  rownames_to_column("cluster") |>
                  pivot_longer(!cluster, names_to = "sample_id", values_to = "CLR_cleanY_Sex"))

if(is_cell_type) {
    ## clr plot of top results
    clr_boxplot_Sex_Micro <- clr_prop_long |>
        filter(grepl("Micro", cluster)) |>
        ggplot(aes(x = Sex, y = CLR_cleanY_Sex, fill = Sex)) +
        geom_boxplot(outlier.shape = NA) +
        # geom_jitter(aes(color = error), width = .1) +
        geom_jitter(width = .1) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = sex_colors) +
        theme_bw() +
        theme(legend.position = "None")

    ggsave(clr_boxplot_Sex_Micro , filename = here(plot_dir, "xenium_clr_boxplot_Sex_Micro.png"), height = 4, width = 6)

    clr_boxplot_Sex_Oligo <- clr_prop_long |>
        filter(grepl("Oligo", cluster)) |>
        ggplot(aes(x = Sex, y = CLR_cleanY_Sex, fill = Sex)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = sex_colors) +
        theme_bw() +
        theme(legend.position = "None")

    ggsave(clr_boxplot_Sex_Oligo, filename = here(plot_dir, "xenium_clr_boxplot_Sex_Oligo.png"), height = 4, width = 6)
}


#### Plot cleanY on cobj Sex + APOE carrier ####

mod <- model.matrix( ~ Sex + APOE_carrier + Age + Anc_Afr + chip , erc_info)
cleanY_cobj_carrier_sex <- cleaningY(cobj$E , mod, P=3)

clr_prop_long <- clr_prop_long |>
    left_join(cleanY_cobj_sex |>
                  as.data.frame() |>
                  rownames_to_column("cluster") |>
                  pivot_longer(!cluster, names_to = "sample_id", values_to = "CLR_cleanY_carrier_Sex")) |>
    mutate(Sex_APOE = paste(Sex, APOE_carrier))

if(is_cell_type) {
    ## plots
    clr_boxplot_Sex_carrier_Oligo <- clr_prop_long |>
        filter(grepl("Oligo", cluster)) |>
        ggplot(aes(x = Sex_APOE, y = CLR_cleanY_carrier_Sex, fill = Sex)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1, aes(color = APOE_carrier)) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = sex_colors) +
        scale_color_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "None",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave(clr_boxplot_Sex_carrier_Oligo, filename = here(plot_dir, "xenium_clr_boxplot_Sex_carrier_Oligo.png"), height = 4, width = 6)

    ## plot prop
    prop_boxplot_Sex_carrier_Oligo <- clr_prop_long |>
        filter(grepl("Oligo", cluster)) |>
        ggplot(aes(x = Sex_APOE, y = prop, fill = Sex)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1, aes(color = APOE_carrier)) +
        facet_wrap(~cluster, nrow = 1) +
        scale_fill_manual(values = sex_colors) +
        scale_color_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "None",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave(prop_boxplot_Sex_carrier_Oligo, filename = here(plot_dir, "xenium_prop_boxplot_Sex_carrier_Oligo.png"), height = 4, width = 6)


    clr_boxplot_Sex_carrier_Astro <- clr_prop_long |>
        filter(grepl("Astro", cluster)) |>
        ggplot(aes(x = Sex_APOE, y = CLR_cleanY_carrier_Sex, fill = Sex)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1, aes(color = APOE_carrier)) +
        facet_wrap(~cluster, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = sex_colors) +
        scale_color_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "None",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave(clr_boxplot_Sex_carrier_Astro, filename = here(plot_dir, "xenium_clr_boxplot_Sex_carrier_Astro.png"), height = 4, width = 6)

    clr_boxplot_Sex_carrier_Micro <- clr_prop_long |>
        filter(grepl("Micro", cluster)) |>
        ggplot(aes(x = Sex_APOE, y = CLR_cleanY_carrier_Sex, fill = Sex)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = .1, aes(color = APOE_carrier)) +
        facet_wrap(~cluster, nrow = 1, scales = "free_y") +
        scale_fill_manual(values = sex_colors) +
        scale_color_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "None",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave(clr_boxplot_Sex_carrier_Micro, filename = here(plot_dir, "xenium_clr_boxplot_Sex_carrier_Micro.png"), height = 4, width = 6)
}

if(is_cell_type) {
    #### compare to snRNA-seq crumblr Sex tests ####
    sn_diff_prop_Sex <- read_csv(here("processed-data", "04_snRNA-seq", "22_crumblr_sn", "sn_diff_prop_tree_test_Sex.csv"))

    diff_prop_tree_Sex_compare <- res_tb |>
        rename_at(5:12, x_name) |>
        select(-branch.length) |> ## stored branch lengths cause buggy join
        full_join(sn_diff_prop_Sex) |>
        mutate(
            text = "Sex",
            label_multi = label |> str_split("/"),
            n_labels = lengths(label_multi),
            label_simple = map_chr(label_multi, ~ .x |>
                                       str_remove_all("\\.[A-Za-z0-9_]+") |>
                                       unique() |>
                                       paste(collapse = "/")
            ),
            label_simple2 = ifelse(n_labels < 3, label, paste0(label_simple, '[', node,']')),
            cell_type_broad = ifelse(grepl("/", label_simple), "MULTI", label_simple),
            signif_cat = case_when(
                pvalue < p_cutoff  & x_pvalue < p_cutoff ~ "both",
                pvalue < p_cutoff  & x_pvalue >= p_cutoff ~ "sn only",
                pvalue >= p_cutoff & x_pvalue < p_cutoff ~ "xenium only",
                pvalue >= p_cutoff & x_pvalue >= p_cutoff ~ "neither"),
            type = ifelse(n_labels > 1, "node", "leaf")
        )

    signif_colors <- c(both = "purple", `xenium only` = "blue", `sn only` = "red")

    tree_compare_beta_scatter <- diff_prop_tree_APOE_carrier_compare |>
        ggplot(aes(beta, x_beta, color = signif_cat)) +
        geom_point() +
        geom_abline(color = "red", linetype = "dashed") +
        ggrepel::geom_text_repel(aes(label = label_simple2), size = 1.7) +
        scale_color_manual(values = signif_colors) +
        theme_bw() +
        labs(x = "snRNA-seq beta", y = "Xenium beta", title = "Crumblr: ~Sex")

    ggsave(tree_compare_beta_scatter,
           filename = here(plot_dir, "tree_compare_beta_scatter_Sex.png"), width = 6, height = 5)

    ggsave(tree_compare_beta_scatter + facet_wrap(~type, nrow = 1),
           filename = here(plot_dir, "tree_compare_beta_scatter_Sex_facet.png"), width = 11, height = 5)
}


#### fit Age ####
message(Sys.time(), " - Fit Age")

# Extract results for each cluster
(diff_prop_Age<- topTable(fit, coef = "Age", number = Inf))

#                           logFC    AveExpr           t    P.Value adj.P.Val         B
# Excit.L6_CT       0.0315704435 -0.6400153  1.97599387 0.07654914 0.9377408 -5.107210
# Astro.4           0.0437897522  1.0228549  1.92446545 0.08336477 0.9377408 -4.762592
# Micro.4          -0.0478180550 -3.3833630 -1.86252019 0.09230820 0.9377408 -5.250839

write.csv(diff_prop_Age, file = here(data_dir, "xenium_diff_prop_Age.csv"))

# Perform multivariate test across the hierarchy
res <- treeTest(fit, cobj, tree.clusCollapsed, coef = "Age")

res_tb <- res |> as_tibble()

res_tb |> arrange(FDR) |> select(label, beta,pvalue, FDR)
res_tb |> write_csv(here(data_dir, "xenium_diff_prop_tree_test_Age.csv"))

# res_tb <- read_csv(here(data_dir, "xenium_diff_prop_tree_test_Age.csv"))

# Plot hierarchy and testing results
nominal_fdr_cut <- res_tb |> filter(pvalue < 0.1) |> slice_max(FDR) |> pull(FDR) |> max()

tree_fdr_plot <- plotTreeTest(res, fdr.cutoff = nominal_fdr_cut)
ggsave(tree_fdr_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_Age_Tree_FDR.png"))

# Plot hierarchy and regression coefficients
tree_beta_plot <- plotTreeTestBeta(res, fdr.cutoff = nominal_fdr_cut) +
    theme(legend.position = "left" 
          # legend.box = "vertical"
    )
ggsave(tree_beta_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_Age_Tree_beta.png"), width = 4, height = 7)

# forest plots
forest_plot <- plotForest(res, hide = FALSE) 
ggsave(forest_plot, filename = here(plot_dir, "xenium_crumblr_cell_type_Age_forest.png"), width = 4, height = 7)

#### Plot cleanY CLR for Age ####
message(Sys.time(), " - cleanY Age")

mod <- model.matrix( ~ Age + APOE_carrier + Anc_Afr + chip , erc_info)
cleanY_cobj_Age <- cleaningY(cobj$E , mod, P=2)

clr_prop_long$CLR_cleanY_Age <- NULL

clr_prop_long <- clr_prop_long |> 
    left_join(cleanY_cobj_Age |>
                  as.data.frame() |>
                  rownames_to_column("cluster") |>
                  pivot_longer(!cluster, names_to = "sample_id", values_to = "CLR_cleanY_Age")) 


#### Plot cleanY CLR for Age & APOE ####
message(Sys.time(), " - plot CLR for AGE")

# mod <- model.matrix( ~ Age + APOE_carrier + Anc_Afr + chip , erc_info)
cleanY_cobj_Age <- jaffelab::cleaningY(cobj$E , mod, P=3)

clr_prop_long$CLR_cleanY_carrier_Age <- NULL

clr_prop_long <- clr_prop_long |> 
    left_join(cleanY_cobj_Age |>
                  as.data.frame() |>
                  rownames_to_column("cluster") |>
                  pivot_longer(!cluster, names_to = "sample_id", values_to = "CLR_cleanY_carrier_Age")) 

if(is_cell_type) {
    clr_sactter_Age_Oligo <- clr_prop_long |>
        filter(grepl("Oligo", cluster)) |>
        ggplot(aes(x = Age, y = CLR_cleanY_Age)) +
        geom_point() +
        geom_smooth(method = "lm") +
        facet_wrap(~cluster, nrow = 1) +
        theme_bw() 

    ggsave(clr_sactter_Age_Oligo, filename = here(plot_dir, "xenium_clr_sactter_Age_Oligo.png"), height = 4, width =10)

    clr_sactter_Age_Oligo_APOE <- clr_prop_long |>
        filter(grepl("Oligo", cluster)) |>
        ggplot(aes(x = Age, y = CLR_cleanY_carrier_Age, color = APOE_carrier)) +
        geom_point() +
        geom_smooth(method = "lm") +
        facet_wrap(~cluster, nrow = 1) +
        scale_color_manual(values = APOE_carrier_colors) +
        theme_bw() +
        theme(legend.position = "bottom")

    ggsave(clr_sactter_Age_Oligo_APOE, filename = here(plot_dir, "xenium_clr_sactter_Age_Oligo_APOE.png"), height = 4, width =10)

    #### compare to snRNA-seq crumblr Age tests ####
    message(Sys.time(), " - compare to snRNA-seq crumblr Age")

    sn_diff_prop_Age <- read_csv(here("processed-data", "04_snRNA-seq", "22_crumblr_sn", "sn_diff_prop_tree_test_Age.csv")) 

    diff_prop_tree_Age_compare <- res_tb |> 
        rename_at(5:12, x_name) |>
        select(-branch.length) |> ## stored branch lengths cause buggy join
        full_join(sn_diff_prop_Age) |>
        mutate(
            test = "Age",
            label_multi = label |> str_split("/"),
            n_labels = lengths(label_multi),
            label_simple = map_chr(label_multi, ~ .x |>
                                       str_remove_all("\\.[A-Za-z0-9_]+") |>
                                       unique() |>
                                       paste(collapse = "/")
            ),
            label_simple2 = ifelse(n_labels < 3, label, paste0(label_simple, '[', node,']')),
            cell_type_broad = ifelse(grepl("/", label_simple), "MULTI", label_simple),
            signif_cat = case_when(
                pvalue < p_cutoff  & x_pvalue < p_cutoff ~ "both",
                pvalue < p_cutoff  & x_pvalue >= p_cutoff ~ "sn only",
                pvalue >= p_cutoff & x_pvalue < p_cutoff ~ "xenium only",
                pvalue >= p_cutoff & x_pvalue >= p_cutoff ~ "neither"),
            type = ifelse(n_labels > 1, "node", "leaf")
        )

    signif_colors <- c(both = "purple", `xenium only` = "blue", `sn only` = "red")

    tree_compare_beta_scatter <- diff_prop_tree_APOE_carrier_compare |>
        ggplot(aes(beta, x_beta, color = signif_cat)) +
        geom_point() +
        geom_abline(color = "red", linetype = "dashed") +
        ggrepel::geom_text_repel(aes(label = label_simple2), size = 1.7) +
        scale_color_manual(values = signif_colors) +
        theme_bw() +
        labs(x = "snRNA-seq beta", y = "Xenium beta", title = "Crumblr: ~Age") 

    ggsave(tree_compare_beta_scatter, 
           filename = here(plot_dir, "tree_compare_beta_scatter_Age.png"), width = 6, height = 5)

    ggsave(tree_compare_beta_scatter + facet_wrap(~type, nrow = 1), 
           filename = here(plot_dir, "tree_compare_beta_scatter_Age_facet.png"), width = 11, height = 5)

    #### export diff_prop_tree_*_compare tables ####

    diff_prop_tree_compare <- rbind(
        diff_prop_tree_APOE_carrier_compare,
        diff_prop_tree_Sex_compare,
        diff_prop_tree_Age_compare
    )

    diff_prop_tree_compare |> write_csv(here(data_dir, "ERC_xenium_sn_diff_prop_tree_compare.csv"))
}


# slurmjobs::job_single('11_xenium_crumblr', create_shell = TRUE, memory = '10G', command = "Rscript 11_xenium_crumblr.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
