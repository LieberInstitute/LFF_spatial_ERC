## Louise Huuki-Myers, Aug 2025
## Find overlap in DEG and LR data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("data.table")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"

data_dir <- here("processed-data", "17_cell_cell_comm", "10_LR_DEG_check")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "17_cell_cell_comm", "10_LR_DEG_check")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)

#### load DE data ####
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn) |> select(gene_name, cluster, vlmf_t, vlmf_adj.P.Val, vlmf_logFC)

DE_data_signif <- DE_data |> filter(vlmf_adj.P.Val < 0.05)
DE_data_signif |> count(cluster)


DE_data |> filter(gene_name == "NRXN1") |> arrange(vlmf_adj.P.Val)
# gene_name cluster   vlmf_t vlmf_adj.P.Val vlmf_logFC
#     1 NRXN1     Oligo.3     3.27         0.0398      1.15 

DE_data |> filter(gene_name == "NLGN1") |> arrange(vlmf_adj.P.Val)
# gene_name cluster      vlmf_t vlmf_adj.P.Val vlmf_logFC
#1 NLGN1     Oligo.3      -2.74          0.0739     -0.538


#top pairs 
LR_top_pairs <- tibble(LR = c('NPTX1^NPTXR', 'RTN4^LINGO1', 'GNAS^ADCY1', 'RTN4^RTN4R', 'CALM1^PTPRA', 'PSAP^SORT1')) |>
    separate(LR, sep= "\\^", into = c("L", "R"), remove = FALSE)

LR_top_pairs_long <- LR_top_pairs |> 
    pivot_longer(!LR, names_to = "class", values_to = "gene_name")
    
LR_top_pairs_long |> left_join(DE_data_signif)
# LR          class gene_name cluster vlmf_t vlmf_adj.P.Val vlmf_logFC
# <chr>       <chr> <chr>     <fct>    <dbl>          <dbl>      <dbl>
# 1 NPTX1^NPTXR L     NPTX1     NA       NA           NA          NA    
# 2 NPTX1^NPTXR R     NPTXR     Astro.1   5.73         0.0310      0.751
# 3 NPTX1^NPTXR R     NPTXR     Oligo.3   4.16         0.0178      1.21 
# 4 RTN4^LINGO1 L     RTN4      NA       NA           NA          NA    
# 5 RTN4^LINGO1 R     LINGO1    NA       NA           NA          NA    
# 6 GNAS^ADCY1  L     GNAS      NA       NA           NA          NA    
# 7 GNAS^ADCY1  R     ADCY1     NA       NA           NA          NA    
# 8 RTN4^RTN4R  L     RTN4      NA       NA           NA          NA    
# 9 RTN4^RTN4R  R     RTN4R     Oligo.3   3.18         0.0443      1.48 
# 10 CALM1^PTPRA L     CALM1     NA       NA           NA          NA    
# 11 CALM1^PTPRA R     PTPRA     NA       NA           NA          NA    
# 12 PSAP^SORT1  L     PSAP      NA       NA           NA          NA    
# 13 PSAP^SORT1  R     SORT1     NA       NA           NA          NA   

#### bivariate data ####
bivariate_data <- fread(here("processed-data","17_cell_cell_comm","liana","bivariate_stats.csv.gz"))

bivariate_data |> filter(ligand == "ADGRB1", receptor == "RTN4R")

bivariate_data_summary <- bivariate_data |> 
    group_by(ligand, receptor) |> 
    summarise(mean_morans = mean(morans),
              median_morans = median(morans),
              mean_mean = mean(mean),
              max_moran_pval = max(morans_pvals),
              n_moran_pval_pass = sum(morans_pvals < 0.05)) |>
    mutate(LR = paste0(ligand, "->", receptor))

bivariate_data_summary |> arrange(-mean_morans)

bivariate_data_summary |> 
    pivot_longer(!c("ligand", "receptor", "LR"), names_to = "metric") |>
    ggplot(aes(x = value)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(~metric, ncol = 1, scales = "free_x")

bivariate_data_summary_scatter <- bivariate_data_summary |>
    ggplot(aes(x = mean_morans, y = mean_mean, color = n_moran_pval_pass)) +
    geom_point() +
    geom_text_repel(aes(label = LR))

bivariate_data_summary_pass_scatter <- bivariate_data_summary |>
    ggplot(aes(x = mean_morans, y = n_moran_pval_pass, color = n_moran_pval_pass)) +
    geom_point(aes())

bivariate_data_summary |>
    ggplot(aes(x = median_morans, y = -log10(max_moran_pval), color = n_moran_pval_pass)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05))

# Moran’s R values near zero imply spatial independence, while positive or negative values reflect spatial co-clustering or spatial cross-dispersion
bivariate_data_summary |> ungroup() |> count(n_moran_pval_pass) |> arrange(-n_moran_pval_pass)

## moran p-val pass all 
bivariate_data_summary |> arrange(-mean_morans)

# Mean "global bivariate score"
bivariate_data_summary |> arrange(-mean_mean)


# top-expressed and highest-spatial-association LR pairs

#### load LR data ####
#pre MFA non-spatially aware data
liana_fn <- here("processed-data", "17_cell_cell_comm", "liana", "ranked_results", sprintf("%s.csv.gz", tolower(gsub("sn_", "", opt$datatype))))
file.exists(liana_fn)

liana_data <- fread(liana_fn)

# rank_aggregate method, which provides a robust rank consensus that combines the predictions of multiple ligand-receptor methods
liana_data |> count(magnitude_rank < 0.05)
# magnitude_rank < 0.05        n
# <lgcl>    <int>
# 1:                 FALSE 12488436
# 2:                  TRUE   783746

liana_data |> count(specificity_rank < 0.05)
# specificity_rank < 0.05        n
# <lgcl>    <int>
# 1:                   FALSE   616979
# 2:                    TRUE    94080
# 3:                      NA 1256112


liana_data |> count(specificity_rank < 0.05, magnitude_rank < 0.05)

summary(liana_data)

liana_data_summary <- liana_data |>
    group_by(source, target, ligand_complex, receptor_complex) |>
    summarise(n_pass_magnitude_rank = sum(magnitude_rank < 0.05),
              n_pass_specificity_rank = sum(specificity_rank < 0.05, na.rm = TRUE),
              n_test = n()) |>
    mutate(source = factor(source, levels = names(cell_type_colors$anno)),
           target = factor(target, levels = names(cell_type_colors$anno)))

summary(liana_data_summary)


liana_data_summary |> ungroup() |> count(n_pass_specificity_rank, n_pass_magnitude_rank) 

liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20) |> count(source, target) |> arrange(-n)
liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20) |> count(ligand_complex, receptor_complex) |> arrange(-n)

liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20, target == "Oligo.3") |> count(source, target) |> arrange(-n)
liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20, source == "Oligo.3") |> count(source, target) |> arrange(-n)

## Oligo.3 <-> astro
liana_data_summary |> filter(source == "Astro.3", target == "Oligo.3") |> arrange(-n_pass_magnitude_rank)
liana_data_summary |> filter(grepl("Astro", source), target == "Oligo.3") |> arrange(-n_pass_magnitude_rank)

liana_data_summary |> filter(n_pass_magnitude_rank >= 20, target == "Oligo.3")


liana_data_summary |> filter(n_pass_magnitude_rank >= 20, ligand_complex == "NLGN1", receptor_complex == "NRXN1")

liana_data_summary |> filter(n_pass_magnitude_rank == 30) |> ungroup() |> count(source)
liana_data_summary |> filter(n_pass_magnitude_rank == 30) |> ungroup() |> count(target)


liana_data_summary |> filter(n_pass_magnitude_rank == 30) 

## plots
liana_data_summary |> ungroup() |> 
    count(n_pass_specificity_rank, n_pass_magnitude_rank)  |>
    filter(n_pass_magnitude_rank > 10,
           n_pass_specificity_rank > 0) |>
    ggplot(aes(x = n_pass_magnitude_rank, n_pass_specificity_rank, fill = n)) +
    geom_tile()

liana_data_summary |> 
    ggplot(aes(x = n_pass_magnitude_rank, y= n_pass_specificity_rank)) +
    geom_point()

n_pass_magnitude_rank_histo <- liana_data_summary |> 
    filter(n_pass_magnitude_rank > 0) |>
    ggplot(aes(x = n_pass_magnitude_rank)) +
    geom_histogram(binwidth = 1) +
    theme_bw()

ggsave(n_pass_magnitude_rank_histo, filename = here(plot_dir, "n_pass_magnitude_rank_histo.png"))

# liana_data_summary |> 
#     filter(n_pass_magnitude_rank > 0) |>
#     ggplot(aes(x = n_pass_magnitude_rank)) +
#     geom_histogram(binwidth = 1) +
#     theme_bw() +
#     facet_wrap(~target)

source_target_counts <- liana_data_summary |> 
    filter(n_pass_magnitude_rank > 20) |>
    ungroup() |>
    count(source, target) |>
    group_by(source) |>
    mutate(p_target = n/sum(n)) |>
    group_by(target) |>
    mutate(p_source = n/sum(n))

LR_count_tile <- source_target_counts |>
    ggplot(aes(x = target, y = source, fill = n)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Number\nLR pairs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(LR_count_tile, filename = here(plot_dir, "LR_count_tile.png"), height = 6, width = 7)

LR_p_source_tile <- source_target_counts |>
    ggplot(aes(x = target, y = source, fill = p_source)) +
    geom_tile() +
    scale_fill_viridis_c(name = "prop source\nLR pairs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(LR_p_source_tile, filename = here(plot_dir, "LR_p_source_tile.png"), height = 6, width = 7)
    
LR_p_target_tile <- source_target_counts |>
    ggplot(aes(x = target, y = source, fill = p_target)) +
    geom_tile() +
    scale_fill_viridis_c(name = "prop target\nLR pairs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
ggsave(LR_p_target_tile, filename = here(plot_dir, "LR_p_target_tile.png"), height = 6, width = 7)


liana_data_summary |> 
    filter(n_pass_magnitude_rank > 20) |>
    ungroup() |>
    count(source, target) |>
    filter(source == "Oligo.3") |> 
    arrange(-n)



