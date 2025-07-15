## Louise Huuki-Myers, July 2025
## Detirmine optimal resolution for DE

#### set-up ####
library("tidyverse")
library("sessioninfo")
library("here")

## set up output dirs
data_dir <- here("processed-data", "12_voomLmFit", "02_optimal_cluster_voomLmFit")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "12_voomLmFit", "02_optimal_cluster_voomLmFit")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)

#### Read summary ####

## broad
broad_summary <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE","sn_broad", "DGE_results_carrier_sn_broad.Rds")) |>
    group_by(cell_type_broad = cluster) |>
    filter(vlmf_P.Value < 0.05) |>
    summarise(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              n_pval05 = sum(vlmf_P.Value < 0.05),
              n_pval05_unique = n_pval05,
              pval05_genes = list(gene_name),
              pval05_genes_unique = pval05_genes) |>
    mutate(hc = "cell_type_broad",
           n= 1)

## fine
fine_count <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE","sn_fine", "DGE_results_carrier_sn_fine.Rds")) |>
    group_by(cell_type_fine = cluster) |>
    filter(vlmf_P.Value < 0.05) |>
    summarise(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              n_pval05 = sum(vlmf_P.Value < 0.05),
              pval05_genes = list(gene_name))  |>
    mutate(cell_type_fine = as.character(cell_type_fine),
           cell_type_broad = jaffelab::ss(cell_type_fine, "\\."),
           hc = "cell_type_fine") 

fine_unique_genes <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE","sn_fine", "DGE_results_carrier_sn_fine.Rds")) |>
    mutate(cell_type_broad = jaffelab::ss(as.character(cluster), "\\.")) |>
    filter(vlmf_P.Value < 0.05) |>
    group_by(cell_type_broad) |>
    summarise(n_pval05_unique = length(unique(gene_name)),
              pval05_genes = list(gene_name),
              pval05_genes_unique = list(unique(gene_name))
              )
    
fine_summary <- fine_count |>
    group_by(cell_type_broad, hc) |>
    summarise(n = n(),
              n_FDR05 = sum(n_FDR05),
              n_pval05 = sum(n_pval05)) |>
    left_join(fine_unique_genes)
    
## hc clusters gene sets
hc_data <- map_dfr(list.files(here(data_dir), recursive = TRUE, pattern = ".rds", full.names = TRUE), ~do.call("rbind", readRDS(.x)))

hc_count <- hc_data |>
    group_by(cell_type_hc = cluster) |>
    filter(P.Value < 0.05) |>
    summarise(n_FDR05 = sum(adj.P.Val < 0.05),
              n_pval05 = sum(P.Value < 0.05),
              pval05_genes = list(gene_name))  |>
    mutate(cell_type_broad = jaffelab::ss(cell_type_hc, "\\."),
           hc = jaffelab::ss(cell_type_hc, "\\.", 2)) 

hc_unique_genes <- hc_data |>
    mutate(cell_type_broad = jaffelab::ss(cluster, "\\."),
           hc = jaffelab::ss(cluster, "\\.", 2))  |>
    filter(P.Value < 0.05) |>
    group_by(cell_type_broad, hc) |>
    summarise(n_pval05_unique = length(unique(gene_name)),
              pval05_genes = list(gene_name),
              pval05_genes_unique = list(unique(gene_name))
    )

hc_summary <- hc_count |>
    group_by(cell_type_broad, hc) |>
    summarise(n = n(),
              n_FDR05 = sum(n_FDR05),
              n_pval05 = sum(n_pval05)) |>
    left_join(hc_unique_genes) |>
    mutate(hc = paste0("cell_type_hc", hc))

## hc clusters

opt_clus_count <- map_dfr(list.files(data_dir, pattern = "optimal_cluster_FDR05_summary", recursive = TRUE, full.names = TRUE), read.csv)

setequal(colnames(hc_summary), colnames(broad_summary))

opt_clus_summary <- hc_summary |>
    bind_rows(broad_summary) |>
    bind_rows(fine_summary) |>
    ungroup() |>
    mutate(hc = fct_relevel(factor(hc), "cell_type_fine", after = Inf))

opt_clus_summary |> count(hc)

sort(unique(opt_clus_summary$hc))

opt_clus_summary |> group_by(cell_type_broad)|> slice_max(n_FDR05) |> select(-pval05_genes, -pval05_genes_unique)
# cell_type_broad hc                  n n_FDR05 n_pval05 n_pval05_unique
# <chr>           <fct>           <dbl>   <int>    <int>           <int>
# 1 Astro           cell_type_fine      5     127     4831            3699
# 2 Excit           cell_type_hc2       1     321     3889            3888
# 3 Inhib           cell_type_hc2       1      34     1950            1950
# 4 Inhib           cell_type_broad     1      34     1950            1950
# 5 Macro           cell_type_broad     1      14      375             375
# 6 Macro           cell_type_fine      1      14      375             375
# 7 Micro           cell_type_hc4       4      68     2574            2096
# 8 OPC             cell_type_hc4       4      10     2303            2005
# 9 Oligo           cell_type_hc2       2    1045     6018            5167
# 10 Vasc            cell_type_fine      2       7     1351            1273

opt_clus_summary |> group_by(cell_type_broad)|> slice_max(n_pval05) 


#### unqiue gene sets ####

## find common pval05 genes
genes_inter <- opt_clus_summary |>
    group_by(cell_type_broad) |>
    summarise(gene_list = list(pval05_genes_unique)) |>
    rowwise()|>
    mutate(pval05_genes_inter = list(Reduce(intersect, gene_list)),
           n_pval05_genes_inter = length(pval05_genes_inter))


## compare sets
opt_clus_summary <- opt_clus_summary |>
    left_join(genes_inter |> select(cell_type_broad, n_pval05_genes_inter)) |>
    mutate(n_pval05_hc_unique = n_pval05_unique - n_pval05_genes_inter)

test2 |>
    ggplot(aes(x = hc, y = n_pval05_hc_unique, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ncell_dot <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n, color = cell_type_broad)) +
    # geom_line(aes(group = cell_type_broad)) +
    geom_jitter(width = 0, height = 0.2) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(ncell_dot, filename = here(plot_dir, "hc_ncell_dot.png"))

## n FDR<05 vs. hc level
nFDR_line <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n_FDR05, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(nFDR_line, filename = here(plot_dir, "hc_nFDR_line.png"))

FDR_pval_scatter <- opt_clus_summary |> 
    ggplot(aes(x = n_pval05, y = n_FDR05, color = cell_type_broad)) +
    # geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(FDR_pval_scatter, filename = here(plot_dir, "hc_FDR_pval_scatter.png"))

## n pval<05 vs. hc level
nPval_line <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n_pval05, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(nPval_line, filename = here(plot_dir, "hc_nPval_line.png"))

## unqiue between same hc groups
nPval_line <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n_pval05_unique, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(nPval_line, filename = here(plot_dir, "hc_nPval_unique_line.png"))

## unique n pval<05 vs. hc level
nPval_hc_unique_line <- opt_clus_summary |>
    ggplot(aes(x = hc, y = n_pval05_hc_unique, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(nPval_hc_unique_line, filename = here(plot_dir, "hc_nPval_hc_unique_line.png"))

## set of genes pval,0.5
nPval_hc_common_col <- genes_inter |>
    ggplot(aes(x = cell_type_broad, y = n_pval05_genes_inter, fill = cell_type_broad)) +
    geom_col() +
    theme_bw() +
    scale_fill_manual(values = cell_type_colors$broad)

ggsave(nPval_hc_common_col, filename = here(plot_dir, "hc_nPval_common_col.png"))


