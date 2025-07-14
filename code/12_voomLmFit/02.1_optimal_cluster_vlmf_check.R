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

broad_summary <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE","sn_broad", "DGE_results_carrier_sn_broad.Rds")) |>
    group_by(cell_type_broad = cluster) |>
    summarise(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
               n_pval05 = sum(vlmf_P.Value < 0.05)) |>
    mutate(hc = "cell_type_broad",
           n= 1)

fine_count <- read.csv(here("processed-data", "12_voomLmFit", "01_Clusterwise_voomLmFit", "vlmf_sn_fine", "vlmf_FDR05_summary-sn_fine.csv")) |>
    dplyr::select(cell_type_fine = cluster, n_FDR05 = carrier) |>
    mutate(cell_type_broad = jaffelab::ss(cell_type_fine, "\\."),
           hc = "cell_type_fine")

fine_count <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE","sn_fine", "DGE_results_carrier_sn_fine.Rds")) |>
    group_by(cell_type_fine = cluster) |>
    summarise(n_FDR05 = sum(vlmf_adj.P.Val < 0.05),
              n_pval05 = sum(vlmf_P.Value < 0.05))  |>
    mutate(cell_type_fine = as.character(cell_type_fine),
           cell_type_broad = jaffelab::ss(cell_type_fine, "\\."),
           hc = "cell_type_fine") 

list.files(here("processed-data", "13_compile_DGE", "01_compile_DGE","sn_broad", "DGE_results_carrier_sn_broad.Rds"))

fine_summary <- fine_count |>
    group_by(cell_type_broad, hc) |>
    summarise(n_FDR05 = sum(n_FDR05),
              n_pval05 = sum(n_pval05),
              n = n())
    

opt_clus_count <- map_dfr(list.files(data_dir, pattern = "optimal_cluster_FDR05_summary", recursive = TRUE, full.names = TRUE), read.csv)

opt_clus_summary <- opt_clus_count |>
    group_by(cell_type_broad, hc) |>
    summarise(n_FDR05 = sum(n_FDR05),
              n_pval05 = sum(n_pval05),
              n = n()) |>
    bind_rows(broad_summary) |>
    bind_rows(fine_summary) |>
    ungroup() |>
    mutate(hc = fct_relevel(factor(hc), "cell_type_fine", after = Inf))

opt_clus_summary |> count(hc)

sort(unique(opt_clus_summary$hc))

opt_clus_summary |> group_by(cell_type_broad)|> slice_max(n_FDR05) 
opt_clus_summary |> group_by(cell_type_broad)|> slice_max(n_pval05) 

# cell_type_broad hc              n_FDR05     n
# <chr>           <fct>             <int> <dbl>
# 1 Astro           cell_type_fine      127     5
# 2 Excit           cell_type_hc2       321     1

# 3 Inhib           cell_type_hc2        34     1
# 4 Inhib           cell_type_broad      34     1

# 5 Macro           cell_type_broad      14     1 * identical no split
# 6 Macro           cell_type_fine       14     1

# 7 Micro           cell_type_hc4        68     4
# 8 OPC             cell_type_hc4        10     4
# 9 Oligo           cell_type_hc2      1045     2 

#10 Vasc            cell_type_fine        7     2 * hc failed...

ncell_dot <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n, color = cell_type_broad)) +
    # geom_line(aes(group = cell_type_broad)) +
    geom_jitter(width = 0, height = 0.2) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(ncell_dot, filename = here(plot_dir, "hc_ncell_dot.png"))


nFDR_line <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n_FDR05, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(nFDR_line, filename = here(plot_dir, "hc_nFDR_line.png"))

nPval_line <- opt_clus_summary |> 
    ggplot(aes(x = hc, y = n_pval05, color = cell_type_broad)) +
    geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(nFDR_line, filename = here(plot_dir, "hc_nPval_line.png"))

FDR_pval_scatter <- opt_clus_summary |> 
    ggplot(aes(x = n_pval05, y = n_FDR05, color = cell_type_broad)) +
    # geom_line(aes(group = cell_type_broad)) +
    geom_point(aes(size = n)) +
    theme_bw() +
    scale_color_manual(values = cell_type_colors$broad) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(FDR_pval_scatter, filename = here(plot_dir, "hc_FDR_pval_scatter.png"))




