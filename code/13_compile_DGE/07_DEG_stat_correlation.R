## Louise Huuki-Myers, July 2025
## Creat correlation plots for DEG stats 

#### Set up ####
# library("tidyverse")
library("here")
library("sessioninfo")
library("corrplot")
library("Hmisc")


data_dir <- here("processed-data", "13_compile_DGE", "07_DEG_stat_correlation")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "07_DEG_stat_correlation")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
data_types <- c("sn_broad", "sn_fine", "Visium")
names(data_types) <- data_types

dge_data <- purrr::map(data_types, ~readRDS(here(here("processed-data", "13_compile_DGE", "01_compile_DGE", .x, sprintf("DGE_results_carrier_%s.Rds", .x)))))
purrr::map(dge_data, dim)

#### make wide data ####
t_stat_wide <- purrr::map(dge_data, ~.x |>
                       select(gene_id, gene_id, cluster, vlmf_t) |>
                       pivot_wider(values_from = "vlmf_t", names_from = "cluster") |> 
                       tibble::column_to_rownames("gene_id")
)

map(t_stat_wide, dim)

logFC_wide <- map(dge_data, ~.x |>
                      select(gene_id, gene_id, cluster, vlmf_logFC) |>
                      pivot_wider(values_from = "vlmf_logFC", names_from = "cluster") |> 
                      tibble::column_to_rownames("gene_id")
)
map(logFC_wide, dim)

#### pair-wise correlations ####
## t-stat
corr_t_stat <- map(t_stat_wide, ~rcorr(as.matrix(.x), type = "pearson"))

pdf(here(plot_dir, "DE_t_stat_corr_plot.pdf"))
map(corr_t_stat, ~try(print(corrplot(.x$r,
                           p.mat = .x$P,
                           type = "upper",
                           order = "hclust",
                           tl.col = "black",
                           tl.srt = 45,
                           sig.level = 0.01,
                           insig = "blank"))))
dev.off()

## logFC
corr_logFC <- map(logFC_wide, ~rcorr(as.matrix(.x), type = "pearson"))

pdf(here(plot_dir, "DE_logFC_corr_plot.pdf"))
map(corr_logFC, ~try(print(corrplot(.x$r,
                           p.mat = .x$P,
                           type = "upper",
                           order = "hclust",
                           tl.col = "black",
                           tl.srt = 45,
                           sig.level = 0.01,
                           insig = "blank"))))
dev.off()

#### datatype v data type ####

common_genes <- intersect(rownames(t_stat_wide$sn_broad), rownames(t_stat_wide$Visium))
length(common_genes)

sn_v_visium <- cor(t_stat_wide$sn_broad[common_genes,], 
    t_stat_wide$Visium[common_genes,], 
    use = "complete.obs")

corrplot(sn_v_visium, 
         tl.col = "black",
         order = "hclust")

# sn_v_visium <- rcorr(as.matrix(t_stat_wide$sn_broad[common_genes,]), 
#                      as.matrix(t_stat_wide$Visium[common_genes,]), 
#                      type = "pearson")
# 
# corrplot(sn_v_visium, order = "hclust")
