## Louise Huuki-Myers, April 2025
## Examine covaraites in dataset for DGE
## adapted from https://github.com/LieberInstitute/dlpfc_asd/blob/a250a1d7e20bd754c5f1186aa96ce0752d55e556/code/08_pseudoBulkDGE_s/02_covariate_analysis.R

library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("BayesSpace")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("dplyr")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data")
#if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Set up dirs ####
plot_dir <- here("plots", "08_pseudoBulkDGE_sn", "01_covariate_analysis")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
# message(Sys.time(), " - Load HDF5 sce")
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC"))

data <- readRDS(here("processed-data", "04_snRNA-seq", "17_sn_model_pseudobulk","sce_pseudobulk-cell_type_anno.rds"))

# Extract the relevant columns from the loaded data
plot_data <- as.data.frame(colData(data)) |>
    select(sample_id, age, diagnosis, BayesSpace, nspots, sum_umi, expr_chrM_ratio, pmi)


# Generate nspots boxplots for each BayesSpace domain
plot_name <- paste0("nspots_by_bayesspace_k", k_nice, ".pdf")
pdf(here(plot_dir, plot_name), width = 8, height = 6)
for (domain in unique(plot_data$BayesSpace)) {
    y_max_nspots <- max(plot_data$nspots) + 20
    domain_data <- plot_data |> filter(BayesSpace == domain)
    plot <- ggboxplot(
        domain_data, x = "diagnosis", y = "nspots", 
        color = "diagnosis", palette = c("blue", "red"), 
        add = "jitter", shape = 19, 
        xlab = "Diagnosis", ylab = "Number of Spots"
    ) + 
        geom_text(
            aes(label = sample_id),
            position = position_jitter(width = 0.1, height = 0.2),
            hjust = 0.5, vjust = 1,
            size = 3
        ) + 
        ggtitle(paste("BayesSpace Domain:", domain)) +
        theme_bw() + 
        theme(legend.position = "none") +
        stat_compare_means(aes(group = diagnosis, label = paste0("p = ", after_stat(p.format))),
                           label.y = y_max_nspots, method = "t.test")
    print(plot)
}
dev.off()

## nspots boxplots with all spatial domains in 1 plot
p_values_nspots <- compare_means(nspots ~ diagnosis, data = plot_data, group.by = "BayesSpace", method = "t.test")
p_values_nspots$fdr <- p.adjust(p_values_nspots$p, method = "fdr")
y_max_nspots <- max(plot_data$nspots) + 50
plot <- ggboxplot(
    plot_data, x = "BayesSpace", y = "nspots", 
    color = "diagnosis", palette = c("blue", "red"), 
    add = "jitter", shape = 19, 
    xlab = "BayesSpace Domain", ylab = "Number of Spots"
) + 
    geom_text(aes(label = sample_id, color = diagnosis), position = position_jitter(width = 0.1, height = 0.2), size = 1.75, hjust = 0.5, vjust = 1) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    stat_compare_means(aes(group = diagnosis,label = paste0("p = ", after_stat(p.format))),label.y = y_max_nspots, method = "t.test")+
    geom_text(data = p_values_nspots, aes(x = BayesSpace, y = y_max_nspots - 20, label = paste0("FDR = ", signif(fdr, 3))), size = 3, color = "black")
plot_name <- paste0("nspots_boxplot_k", k_nice, ".png")  
ggsave(filename = here(plot_dir, plot_name), plot = plot, width = 12, height = 8)

# Generate sum_umi boxplots across BayesSpace domains
p_values_sum_umi <- compare_means(sum_umi ~ diagnosis, data = plot_data, group.by = "BayesSpace", method = "t.test")
p_values_sum_umi$fdr <- p.adjust(p_values_sum_umi$p, method = "fdr")
y_max_sum_umi <- max(plot_data$sum_umi) + 600000
plot <- ggboxplot(
    plot_data, x = "BayesSpace", y = "sum_umi", 
    color = "diagnosis", palette = c("blue", "red"), 
    add = "jitter", shape = 19, 
    xlab = "BayesSpace Domain", ylab = "sum_umi"
) + 
    geom_text(aes(label = sample_id, color = diagnosis), 
              position = position_jitter(width = 0.1, height = 0.2), 
              size = 1.75, hjust = 0.5, vjust = 1) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    stat_compare_means(aes(group = diagnosis, label = paste0("p = ", after_stat(p.format))),
                       label.y = y_max_sum_umi, method = "t.test") +
    geom_text(data = p_values_sum_umi, aes(x = BayesSpace, y = y_max_sum_umi - 200000, 
                                           label = paste0("FDR = ", signif(fdr, 3))), 
              size = 3, color = "black")

plot_name <- paste0("sum_umi_boxplot_k", k_nice, ".png")
ggsave(filename = here(plot_dir, plot_name), plot = plot, width = 12, height = 8)

# Generate expr_chrM_ratio boxplots with all spatial domains in 1 plot
p_values_mito <- compare_means(expr_chrM_ratio ~ diagnosis, data = plot_data, group.by = "BayesSpace", method = "t.test")
p_values_mito$fdr <- p.adjust(p_values_mito$p, method = "fdr")
y_max_mito <- max(plot_data$expr_chrM_ratio) + 0.02
plot <- ggboxplot(
    plot_data, x = "BayesSpace", y = "expr_chrM_ratio", 
    color = "diagnosis", palette = c("blue", "red"), 
    add = "jitter", shape = 19, 
    xlab = "BayesSpace Domain", ylab = "expr_chrM_ratio"
) + 
    geom_text(aes(label = sample_id, color = diagnosis), 
              position = position_jitter(width = 0, height = 0), 
              size = 1.75, hjust = 0, vjust = 1) + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    stat_compare_means(aes(group = diagnosis, label = paste0("p = ", after_stat(p.format))),
                       label.y = y_max_mito, method = "t.test") +
    geom_text(data = p_values_mito, aes(x = BayesSpace, y = y_max_mito - 0.005, 
                                        label = paste0("FDR = ", signif(fdr, 3))), 
              size = 3, color = "black")

plot_name <- paste0("expr_chrM_ratio_boxplot_k", k_nice, ".png")
ggsave(filename = here(plot_dir, plot_name), plot = plot, width = 12, height = 8)



###### Additional analysis (since age and PMI remain constant) #############
k_nice <- "12"
data_file <- file.path(data_dir, paste0("summed_k", k_nice, ".rds"))
plot_data <- as.data.frame(colData(data)) |>
    select(sample_id, age, diagnosis, BayesSpace, nspots, sum_umi, expr_chrM_ratio, pmi)

plot_data_summary <- plot_data |>
    select(sample_id, age, diagnosis, expr_chrM_ratio, pmi) |>
    distinct(sample_id, .keep_all = TRUE)

# Function for creating and saving boxplots
create_boxplot <- function(data, x, y, y_label, filename) {
    pdf(here(plot_dir, filename), width = 8, height = 6)
    y_max <- max(data[[y]]) + 5
    #print(y_max)
    print(ggboxplot(
        data, x = x, y = y,
        color = x, palette = c("blue", "red"),
        add = "jitter", shape = 19,
        xlab = "Diagnosis", ylab = y_label
    ) +
        geom_text_repel(
            aes(label = sample_id),
            size = 3, color = 'black',
            max.overlaps = Inf,
            force_pull = 2,
            show.legend = FALSE
        ) +
        theme_bw() +
        theme(legend.position = "bottom")+
        stat_compare_means(aes(group = diagnosis,label = paste0("p = ", after_stat(p.format))),label.y = y_max, method = "t.test")
    )
    dev.off()
}

# Generate boxplots
create_boxplot(plot_data_summary, "diagnosis", "pmi", "PMI", "full_pmi_distribution.pdf")
create_boxplot(plot_data_summary, "diagnosis", "age", "Age", "full_age_distribution_by_diagnosis.pdf")

# Function for Violin Plots
create_violin_plot <- function(data, x, y,x_label, y_label, filename) {
    pdf(here(plot_dir, filename), width = 8, height = 6)
    y_max <- max(data[[y]]) + 20
    print(
        ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = diagnosis)) +
            geom_violin(trim = FALSE, alpha = 0.6, drop = FALSE) +
            geom_jitter(aes(color = diagnosis), width = 0.2, alpha = 0.7, size = 2) +
            scale_fill_manual(values = c("blue", "red")) +
            scale_color_manual(values = c("black", "darkgreen")) +
            labs(x = x_label, y = y_label) +
            theme_bw() +
            theme(legend.position = "none")+
            stat_compare_means(aes(group = diagnosis,label = paste0("p = ", after_stat(p.format))),label.y = y_max, method = "t.test")
    )
    dev.off()
}

# Generate Violin Plots
create_violin_plot(plot_data_summary, "diagnosis", "pmi","Diagnosis", "PMI", "pmi_violin.pdf") #y_max = max(data[[y]]) + 20
create_violin_plot(plot_data_summary, "diagnosis", "age","Diagnosis", "AGE", "age_violin.pdf") #y_max <- max(data[[y]]) + 5

# plot_data_summary <- plot_data_summary %>%
#   mutate(age_group = ifelse(age <="21-29", "30-38"))
# # table(plot_data_summary$age_group)
# # 21-29 30-38
# #     9     5
#
# create_violin_plot(plot_data_summary, "age_group", "expr_chrM_ratio", "Age group", "Mitochondrial Expression Ratio", "age_vs_expr_chrM_violin.pdf")
# create_violin_plot(plot_data_summary, "age_group", "pmi", "Age group", "PMI", "age_vs_pmi_violin.pdf")

# Age vs. expr_chrM_ratio scatter plot
pdf(here(plot_dir, "age_vs_expr_chrM_by_diagnosis.pdf"), width = 8, height = 6)
ggplot(plot_data_summary, aes(x = age, y = expr_chrM_ratio, color = diagnosis)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, aes(fill = diagnosis), alpha = 0.2) +
    geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 15) +
    scale_color_manual(values = c("Control" = "blue", "Autism" = "red")) +
    scale_fill_manual(values = c("Control" = "blue", "Autism" = "red")) +
    labs(x = "Age", y = "expr_chrM_ratio", color = "Diagnosis") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()

# #################exploring the mito genes###############
# 
# reference_gtf <- "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
# reference_gtf <- reference_gtf[file.exists(reference_gtf)]
# 
# gtf <- rtracklayer::import(reference_gtf)
# gtf <- gtf[gtf$type == "gene"]
# names(gtf) <- gtf$gene_id
# 
# # Extract mitochondrial genes from GTF file
# mito_genes_gtf <- gtf[seqnames(gtf) == "chrM"]
# 
# # Get the gene IDs and gene names of mitochondrial genes from GTF
# mito_gene_info <- data.frame(
#   gene_id = names(mito_genes_gtf),
#   gene_name = mcols(mito_genes_gtf)$gene_name,
#   seqnames = seqnames(mito_genes_gtf)
# )
# 
# # Find which genes in the spe object correspond to mitochondrial genes
# mito_genes_spe <- rowData(spe)$gene_id %in% mito_gene_info$gene_id
# 
# ##########counts data###############################
# # Extract the counts assay
# counts_data <- assays(spe)$counts
# # Mean counts for each mitochondrial gene
# mean_counts <- rowMeans(counts_data[mito_genes_spe, , drop = FALSE])
# # Sum counts for each mitochondrial gene
# sum_counts <- rowSums(counts_data[mito_genes_spe, , drop = FALSE])
# 
# # Create a dataframe of mitochondrial gene information in the spe object
# mito_genes_spe_info <- data.frame(
#     gene_id = rowData(spe)$gene_id[mito_genes_spe],
#     gene_name = rowData(spe)$gene_name[mito_genes_spe],
#     mean_counts = mean_counts, 
#     sum_counts = sum_counts,
#     seqnames = seqnames(spe)[mito_genes_spe]
# )
# 
# # Reset row names to numeric indexing
# rownames(mito_genes_spe_info) <- NULL
# 
# # Count mitochondrial genes
# mito_gene_count <- nrow(mito_genes_spe_info)
# 
# mito_gene_count
# mito_genes_spe_info
# 
# summary(mito_genes_spe_info)
# 
# 
# # Find which genes in the summed_k09 object correspond to mitochondrial genes
# mito_genes_summed_k09 <- rowData(summed_k09)$gene_id %in% mito_gene_info$gene_id
# 
# # Extract the counts assay
# counts_data_summed_k09 <- assays(summed_k09)$counts
# # Mean counts for each mitochondrial gene
# mean_counts_summed_k09 <- rowMeans(counts_data_summed_k09[mito_genes_summed_k09, , drop = FALSE])
# # Sum counts for each mitochondrial gene
# sum_counts_summed_k09 <- rowSums(counts_data_summed_k09[mito_genes_summed_k09, , drop = FALSE])
# 
# # Create a dataframe of mitochondrial gene information in the spe object
# mito_genes_summed_k09_info <- data.frame(
#     gene_id = rowData(summed_k09)$gene_id[mito_genes_summed_k09],
#     gene_name = rowData(summed_k09)$gene_name[mito_genes_summed_k09],
#     mean_counts = mean_counts_summed_k09, 
#     sum_counts = sum_counts_summed_k09, 
#     seqnames = seqnames(summed_k09)[mito_genes_summed_k09]
# )
# 
# # Reset row names to numeric indexing
# rownames(mito_genes_summed_k09_info) <- NULL
# 
# # Count mitochondrial genes
# mito_gene_count_summed_k09 <- nrow(mito_genes_summed_k09_info)
# 
# mito_gene_count_summed_k09
# mito_genes_summed_k09_info
# 
# summary(mito_genes_summed_k09_info)
# 
# ############logcounts######################################
# 
# # Extract the counts assay
# logcounts_data_summed_k09 <- assays(summed_k09)$logcounts
# # Mean counts for each mitochondrial gene
# mean_logcounts_summed_k09 <- rowMeans(logcounts_data_summed_k09[mito_genes_summed_k09, , drop = FALSE])
# # Sum counts for each mitochondrial gene
# sum_logcounts_summed_k09 <- rowSums(logcounts_data_summed_k09[mito_genes_summed_k09, , drop = FALSE])
# 
# # Create a dataframe of mitochondrial gene information in the summed_test object
# mito_genes_summed_k09_info <- data.frame(
#     gene_id = rowData(summed_k09)$gene_id[mito_genes_summed_k09],
#     gene_name = rowData(summed_k09)$gene_name[mito_genes_summed_k09],
#     mean_logcounts = mean_logcounts_summed_k09, 
#     sum_logcounts = sum_logcounts_summed_k09, 
#     seqnames = seqnames(summed_k09)[mito_genes_summed_k09]
# )
# 
# # Reset row names to numeric indexing
# rownames(mito_genes_summed_k09_info) <- NULL
# 
# # Count mitochondrial genes
# mito_gene_count_summed_test <- nrow(mito_genes_summed_k09_info)
# 
# mito_gene_count_summed_k09
# mito_genes_summed_k09_info
# 
# summary(mito_genes_summed_k09_info)
# 
# # Extract the logcounts assay from spe
# logcounts_data <- assays(spe)$logcounts
# # Mean counts for each mitochondrial gene
# mean_logcounts <- rowMeans(logcounts_data[mito_genes_spe, , drop = FALSE])
# # Sum counts for each mitochondrial gene
# sum_logcounts <- rowSums(logcounts_data[mito_genes_spe, , drop = FALSE])
# 
# # Create a dataframe of mitochondrial gene information in the spe object
# mito_genes_spe_info <- data.frame(
#     gene_id = rowData(spe)$gene_id[mito_genes_spe],
#     gene_name = rowData(spe)$gene_name[mito_genes_spe],
#     mean_logcounts = mean_logcounts, 
#     sum_logcounts = sum_logcounts,
#     seqnames = seqnames(spe)[mito_genes_spe]
# )
# 
# # Reset row names to numeric indexing
# rownames(mito_genes_spe_info) <- NULL
# 
# # Count mitochondrial genes
# mito_gene_count <- nrow(mito_genes_spe_info)
# 
# mito_gene_count
# mito_genes_spe_info
# 
# summary(mito_genes_spe_info)

# slurmjobs::job_single('01_create_pseudobulk_data', create_shell = TRUE, memory = '25G', command = "Rscript 01_create_pseudobulk_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()