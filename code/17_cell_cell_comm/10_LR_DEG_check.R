## Louise Huuki-Myers, Aug 2025
## Find overlap in DEG and LR data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("data.table")
library("ComplexHeatmap")
library("ggrepel")

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

load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

if(opt$datatype == "sn_broad"){
    cluster_colors <- cell_type_colors$broad
    cluster_levels <- names(cell_type_colors$broad)
}else if(opt$datatype == "sn_fine"){
    cluster_colors <- cell_type_colors$anno
    cluster_levels <- names(cell_type_colors$anno)
    
    broad_cell_types <- names(cell_type_colors$broad)
    broad_cell_types <- broad_cell_types[broad_cell_types != "Other"]
    
}else if(opt$datatype == "Visium"){

    cluster_colors <- SpD_colors
    cluster_levels <- names(SpD_colors)
}

cluster_levels <- cluster_levels[cluster_levels != "Other"]

#### load DE data ####
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)

DE_data <- readRDS(DE_data_fn) |> select(gene_name, cluster, vlmf_t, vlmf_adj.P.Val, vlmf_logFC)

DE_data_signif <- DE_data |> filter(vlmf_adj.P.Val < 0.05)
DE_data_signif |> count(cluster)


DE_data |> filter(grepl("NPTX", gene_name), vlmf_adj.P.Val < 0.05) |> arrange(vlmf_adj.P.Val)

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

#### sample level bivariate data ####
bivariate_data <- fread(here("processed-data","17_cell_cell_comm","liana","bivariate_stats.csv.gz"))

bivariate_data |> filter(ligand == "ADGRB1", receptor == "RTN4R")

bivariate_data_summary <- bivariate_data |> 
    group_by(ligand, receptor) |> 
    summarise(mean_morans = mean(morans),
              median_morans = median(morans),
              mean_mean = mean(mean),
              max_moran_pval = max(morans_pvals),
              n_moran_pval_pass = sum(morans_pvals < 0.05)) |>
    ungroup() |>
    arrange(-mean_mean) |>
    mutate(LR = paste0(ligand, "->", receptor),
           mean_rank = row_number())

bivariate_data_summary |> arrange(-mean_morans)
bivariate_data_summary |> arrange(-mean_mean)

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

bivariate_data_summary_pass_scatter_log10 <-bivariate_data_summary |>
    ggplot(aes(x = median_morans, y = -log10(max_moran_pval), color = n_moran_pval_pass)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05))


bivariate_data_summary |> 
    ungroup() |>
    arrange(-mean_mean) |>
    mutate(rank = row_number())

# Moran’s R values near zero imply spatial independence, while positive or negative values reflect spatial co-clustering or spatial cross-dispersion
bivariate_data_summary |> ungroup() |> count(n_moran_pval_pass) |> arrange(-n_moran_pval_pass)

## moran p-val pass all 
bivariate_data_summary |> arrange(-mean_morans)

# Mean "global bivariate score"
bivariate_data_summary |> arrange(-mean_mean)


# top-expressed and highest-spatial-association LR pairs


#### Spot level bivaritate data ####
bivariate_spot_data <- fread(here("processed-data", "17_cell_cell_comm", "liana", "bivariate_stats_spot.csv.gz"))

dim(bivariate_spot_data)

bivariate_spot_data |> 
    group_by(pair_id) |>
    count()

bivariate_SpD_data <- bivariate_spot_data |>
    group_by(pair_id, SpD) |>
    summarise(mean_score = mean(local_score),
              n = n()) |>
    mutate(SpD = factor(SpD, levels = names(SpD_colors)))


bivaraite_SpD_score_histo <- bivariate_SpD_data |>
    ggplot(aes(x = mean_score)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(~SpD) +
    theme_bw()

ggsave(bivaraite_SpD_score_histo, filename = here(plot_dir, "bivariate_SpD_score_histogram.png"))


#### load & Summarize LR data ####
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

## summarize liana data
liana_data_summary <- liana_data |>
    group_by(source, target, ligand_complex, receptor_complex) |>
    summarise(n_pass_magnitude_rank = sum(magnitude_rank < 0.05),
              n_pass_specificity_rank = sum(specificity_rank < 0.05, na.rm = TRUE),
              n_test = n()) |>
    mutate(source = factor(source, levels = cluster_levels),
           target = factor(target, levels = cluster_levels)) |>
    left_join(bivariate_data_summary |> select(ligand_complex = ligand, receptor_complex = receptor, mean_mean, mean_rank))

summary(liana_data_summary)


liana_data_summary |> ungroup() |> count(n_pass_specificity_rank, n_pass_magnitude_rank) 

liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20) |> count(source, target) |> arrange(-n)
liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20) |> count(ligand_complex, receptor_complex) |> arrange(-n)

if(opt$datatype == "sn_fine"){
    
    liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20, target == "Oligo.3") |> count(source, target) |> arrange(-n)
    liana_data_summary |> ungroup() |> filter(n_pass_magnitude_rank >= 20, source == "Oligo.3") |> count(source, target) |> arrange(-n)
    
    ## Oligo.3 <-> astro
    liana_data_summary |> filter(source == "Astro.3", target == "Oligo.3") |> arrange(-n_pass_magnitude_rank)
    liana_data_summary |> filter(grepl("Astro", source), target == "Oligo.3") |> arrange(-n_pass_magnitude_rank)
    
    liana_data_summary |> filter(n_pass_magnitude_rank >= 20, target == "Oligo.3")|> ungroup() |> count(source) |> arrange(-n)
    # source               n
    # <fct>            <int>
    # 1 Excit.L2            15
    # 2 Excit.L6b           14
    # 3 Inhib.Lamp5_Lhx6    11
    # 4 Inhib.Pvalb         11
    
    liana_data_summary |> filter(n_pass_magnitude_rank >= 20, source == "Oligo.3")|> ungroup() |> count(target) |> arrange(-n)
    # target               n
    # <fct>            <int>
    # 1 Excit.L2            24
    # 2 Excit.L6b           17
    # 3 Excit.L5.1          15
    # 4 Inhib.Sst            9
    # 5 Excit.L2_5.1         8
    
    
    liana_data_summary |> filter(n_pass_magnitude_rank >= 20) |> ungroup() |> count(target == "Oligo.3", source == "Oligo.3")
    liana_data_summary |> filter(n_pass_magnitude_rank >= 20) |> ungroup() |> count(target == "Oligo.3", source == "Oligo.3")
    liana_data_summary |> ungroup() |> filter(target == "Oligo.3", source == "Oligo.3") |> arrange(-n_pass_magnitude_rank)
    
    }
   
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

n_pass_magnitude_rank_histo <- liana_data_summary |> 
    filter(n_pass_magnitude_rank > 0) |>
    ggplot(aes(x = n_pass_magnitude_rank)) +
    geom_histogram(binwidth = 1) +
    theme_bw()

ggsave(n_pass_magnitude_rank_histo, filename = here(plot_dir, sprintf("n_pass_magnitude_rank_histo_%s.png", opt$datatype)))

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

ggsave(LR_count_tile, filename = here(plot_dir, sprintf("LR_count_tile_%s.png", opt$datatype)), height = 6, width = 7)

LR_p_source_tile <- source_target_counts |>
    ggplot(aes(x = target, y = source, fill = p_source)) +
    geom_tile() +
    scale_fill_viridis_c(name = "prop source\nLR pairs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(LR_p_source_tile, filename = here(plot_dir, sprintf("LR_p_source_tile_%s.png", opt$datatype)), height = 6, width = 7)

LR_p_target_tile <- source_target_counts |>
    ggplot(aes(x = target, y = source, fill = p_target)) +
    geom_tile() +
    scale_fill_viridis_c(name = "prop target\nLR pairs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(LR_p_target_tile, filename = here(plot_dir, sprintf("LR_p_target_tile_%s.png", opt$datatype)), height = 6, width = 7)

## LR mean_mean
LR_bivaraite_rank_magnitude_pass_histo <- liana_data_summary |>
    filter(n_pass_magnitude_rank > 20) |>
    ggplot(aes(x = mean_mean)) +
    geom_histogram(binwidth = 0.02) +
    theme_bw()

ggsave(LR_bivaraite_rank_magnitude_pass_histo, filename = here(plot_dir, sprintf("LR_bivaraite_rank_magnitude_pass_histo_%s.png", opt$datatype)), height = 6, width = 7)

if(opt$datatype == "sn_fine"){
    ## Oligo.3 senders and targets 
    
    source_target_counts |>
        filter(source == "Oligo.3") |> 
        arrange(-n)
    
    source_target_counts |>
        filter(source == "Oligo.3" & target == "Oligo.3")  
    
    source_target_counts |>
        filter(source == "Oligo.3" | target == "Oligo.3") |>
        mutate(st = paste(source, target)) |>
        select(st, source, target, n)
    
    LR_count_bar_Oligo.3 <- source_target_counts |>
        filter(source == "Oligo.3" | target == "Oligo.3") |> 
        mutate(st = paste(source, target)) |>
        select(st, source, target, n) |>
        pivot_longer(!c(n, st), names_to = "class", values_to = "cell_type") |>
        ggplot(aes(x = cell_type, y = n, fill = cell_type)) +
        geom_col() +
        scale_fill_manual(values = cell_type_colors$anno) +
        facet_wrap(~class, ncol = 1) +
        theme_bw() +
        labs(y="n LR pairs with Oligo.3") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "None")
    
    ggsave(LR_count_bar_Oligo.3, filename = here(plot_dir, "LR_count_bar_Oligo.3.png"), height = 4, width = 6)
    
}else if(opt$datatype == "sn_broad"){
    ## Oligo.3 senders and targets 
    
    source_target_counts |>
        filter(source == "Oligo") |> 
        arrange(-n)
    
    LR_count_bar_Oligo <- source_target_counts |>
        filter(source == "Oligo" | target == "Oligo") |>
        mutate(st = paste(source, target)) |>
        select(st, source, target, n) |>
        pivot_longer(!c(n, st), names_to = "class", values_to = "cell_type") |>
        ggplot(aes(x = cell_type, y = n, fill = cell_type)) +
        geom_col() +
        scale_fill_manual(values = cluster_colors) +
        facet_wrap(~class, ncol = 1) +
        theme_bw() +
        labs(y="n LR pairs with Oligo") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "None")
    
    ggsave(LR_count_bar_Oligo, filename = here(plot_dir, "LR_count_bar_Oligo.png"), height = 4, width = 6)
    
}


#### DEG intersect ####

LR_DEGS <- liana_data_summary |> 
    filter(n_pass_magnitude_rank > 20) |>
    left_join(DE_data |> select(source = cluster, ligand_complex = gene_name, ligand_t = vlmf_t, ligand_adj.P.Val = vlmf_adj.P.Val, ligand_logFC = vlmf_logFC))|>
    left_join(DE_data |> select(target = cluster, receptor_complex = gene_name, receptor_t = vlmf_t, receptor_adj.P.Val = vlmf_adj.P.Val, receptor_logFC = vlmf_logFC)) |>
    mutate(L_DEG = ligand_adj.P.Val < 0.05,
           R_DEG = receptor_adj.P.Val < 0.05)

LR_DEGS |> filter(is.na(R_DEG)) |> select(source, target, ligand_complex, receptor_complex, ligand_adj.P.Val, receptor_adj.P.Val, L_DEG, R_DEG)

LR_DEGS |> ungroup() |> count(L_DEG, R_DEG)

if(any(LR_DEGS$L_DEG | LR_DEGS$R_DEG, na.rm = TRUE)){
    
    LR_DEGS |> ungroup() |> filter(L_DEG) |>  count(source)
    LR_DEGS |> ungroup() |> filter(R_DEG) |>  count(target)
    
    LR_DEGS |> 
        ungroup() |> 
        filter(ligand_complex == "NRXN1",
               receptor_complex == "NLGN1") |> 
        # filter(target == "Oligo.3") |>
        select(source, target, ligand_complex, receptor_complex, ligand_adj.P.Val, receptor_adj.P.Val, L_DEG, R_DEG, mean_mean)
    
    ## log FC heatmap
    
    ligand_degs <- LR_DEGS |> ungroup() |> filter(L_DEG) |> pull(ligand_complex) |> unique()
    receptor_of_ligand_degs <- LR_DEGS |> ungroup() |> filter(L_DEG) |> pull(receptor_complex) |> unique()
    
    receptor_degs <- LR_DEGS |> ungroup() |> filter(R_DEG) |> pull(receptor_complex) |> unique()
    
    source(here("code", "13_compile_DGE", "logFC_heatmap.R"))
    
    cluster_levels = names(cell_type_colors$anno)
    
    LR_DEGS_plot_data <- LR_DEGS |>
        filter(L_DEG | R_DEG) |>
        ungroup() |>
        mutate(interaction = paste(ligand_complex, "->", receptor_complex),
               DEG = case_when(L_DEG ~ "L_DEG",
                               R_DEG ~ "R_DEG",
                               TRUE ~ "Other")) |>
        select(DEG, interaction, source, target) |>
        pivot_longer(!c(DEG, interaction), names_to = "class", values_to = "cell_type") |>
        mutate(class = case_when(class == "target" ~ "Receptor",
                                 class == "source" ~ "Ligand",
                                 TRUE ~ NA
        )) |>
        unique()
    
    LR_DEGS_plot_data |> count(DEG)
    
    LR_DEGS_plot_data2 <- LR_DEGS |>
        filter(L_DEG | R_DEG) |>
        ungroup() |>
        mutate(interaction = paste(ligand_complex, "->", receptor_complex),) |>
        select(interaction, mean_mean, receptor_complex, ligand_complex) |>
        pivot_longer(!c(interaction, mean_mean), names_to = "class", values_to = "gene_name") |>
        mutate(class = case_when(class == "receptor_complex" ~ "Receptor",
                                 class == "ligand_complex" ~ "Ligand",
                                 TRUE ~ NA
        )) |>
        unique()
    
    
    LR_DEGS_plot_data3 <- LR_DEGS_plot_data |>
        left_join(LR_DEGS_plot_data2) |>
        left_join(DE_data |> 
                      select(cell_type = cluster, gene_name, vlmf_t, vlmf_adj.P.Val, vlmf_logFC)) |>
        replace_na(list(mean_mean = 0)) |>
        mutate(interaction = fct_reorder(interaction, mean_mean),
               signif = case_when(vlmf_adj.P.Val < 0.005 ~ "***",
                                  vlmf_adj.P.Val < 0.01 ~"**",
                                  vlmf_adj.P.Val < 0.05 ~"*",
                                  TRUE~"")
        )
    
    LR_DEGS_plot_data3 |> count(DEG, interaction)
    
    max_abs = max(abs(LR_DEGS_plot_data3$vlmf_logFC))
    
    
    LR_DEGS_plot <- LR_DEGS_plot_data3 |>
        ggplot(aes(x = cell_type, y = interaction, fill = vlmf_logFC)) +
        facet_grid(DEG~class, scales = "free", space = "free") +
        geom_tile(color = "black") +
        geom_text(aes(label = signif)) +
        scale_fill_gradient2(high = APOE_carrier_colors[["E4+"]], low = APOE_carrier_colors[["E2+"]])  +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(LR_DEGS_plot, filename = here(plot_dir, sprintf("LR_DEGS_%s.png", opt$datatype)), width = 8, height = 6)
    
    
    LR_DEGS_plot_ct <- LR_DEGS_plot_data3 |>
        ggplot(aes(x = cell_type, y = interaction, fill = cell_type)) +
        facet_grid(DEG~class, scales = "free", space = "free") +
        geom_tile(color = "black") +
        geom_text(aes(label = signif)) +
        scale_fill_manual(values = cluster_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
    
    ggsave(LR_DEGS_plot_ct, filename = here(plot_dir, sprintf("LR_DEGS_ct_%s.png", opt$datatype)), width = 8, height = 6)
    
    
    LR_DEGS_plot_data3 |>
        group_by(interaction) |>
        slice(1)|>
        arrange(mean_mean)
    
    
    LR_DEGS_mean_mean <- LR_DEGS_plot_data3 |>
        group_by(interaction) |>
        slice(1)|>
        ggplot(aes(x = mean_mean, y = interaction)) +
        geom_col() +
        theme_bw()
    
    ggsave(LR_DEGS_mean_mean, filename = here(plot_dir, sprintf("LR_DEGS_mean_mean_%s.png", opt$datatype)), width = 3, height = 6)
    
}

#### Bivaraite DEG data ####

LR_DEGS_SpD_bar <- bivariate_SpD_data |>
    mutate(interaction = gsub("\\^", " -> ", pair_id)) |>
    filter(interaction %in% LR_DEGS_plot_data3$interaction) |>
    ggplot(aes(y = interaction, x = mean_score, fill = SpD)) +
    geom_col(position = "dodge") +
    theme_bw() +
    scale_fill_manual(values = SpD_colors)

ggsave(LR_DEGS_SpD_bar, filename = here(plot_dir, "LR_DEGS_bivarite_SpD_bar.png"))


LR_DEGS_list <- LR_DEGS_plot_data |>
    select(interaction, DEG) |>
    unique() 

LR_DEGS_SpD_tile <- bivariate_SpD_data |>
    mutate(interaction = gsub("\\^", " -> ", pair_id)) |>
    right_join(LR_DEGS_list) |>
    filter(!is.na(mean_score)) |>
    mutate(interaction = factor(interaction, levels = levels(LR_DEGS_plot_data3$interaction))) |>
    ggplot(aes(y = interaction, x = SpD, fill = mean_score)) +
    geom_tile() +
    theme_bw() +
    facet_grid(DEG~., scales = "free_y", space = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(LR_DEGS_SpD_tile, filename = here(plot_dir, "LR_DEGS_bivarite_SpD_tile.png"))


#### AD Risk Genes ####

AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 

## no commonly observed APOE LS-RTs
liana_data_summary |> 
    filter((ligand_complex  == "APOE" | receptor_complex == "APOE")) |>
    arrange(-n_pass_magnitude_rank)

LR_risk <- liana_data_summary |> 
    filter(n_pass_magnitude_rank > 20 & (ligand_complex %in% AD_risk$symbol | receptor_complex %in% AD_risk$symbol)) |>
    left_join(DE_data |> select(source = cluster, ligand_complex = gene_name, ligand_t = vlmf_t, ligand_adj.P.Val = vlmf_adj.P.Val, ligand_logFC = vlmf_logFC))|>
    left_join(DE_data |> select(target = cluster, receptor_complex = gene_name, receptor_t = vlmf_t, receptor_adj.P.Val = vlmf_adj.P.Val, receptor_logFC = vlmf_logFC)) |>
    mutate(L_DEG = ligand_adj.P.Val < 0.05,
           R_DEG = receptor_adj.P.Val < 0.05,
           L_risk = ligand_complex %in% AD_risk$symbol ,
           R_risk = receptor_complex %in% AD_risk$symbol,
           interaction = paste(ligand_complex, "->", receptor_complex)) |>
    ungroup()

# Check for risk genes in commonly observed LR interactions
LR_risk |> count(L_DEG, R_DEG)

## no LR both risk pairs
LR_risk |> count(L_risk, R_risk)

LR_risk_plot_data <- LR_risk |>
    mutate(Risk = case_when(L_risk ~ "L_risk",
                           R_risk ~ "R_risk",
                           TRUE ~ "Other")) |>
    select(Risk, interaction, source, target) |>
    pivot_longer(!c(Risk, interaction), names_to = "class", values_to = "cell_type") |>
    mutate(class = case_when(class == "target" ~ "Receptor",
                             class == "source" ~ "Ligand",
                             TRUE ~ NA
    )) |>
    unique()

LR_risk_plot_data |> count(Risk, cell_type)

LR_risk_plot_data2 <- LR_risk |>
    select(interaction, mean_mean, receptor_complex, ligand_complex) |>
    pivot_longer(!c(interaction, mean_mean), names_to = "class", values_to = "gene_name") |>
    mutate(class = case_when(class == "receptor_complex" ~ "Receptor",
                             class == "ligand_complex" ~ "Ligand",
                             TRUE ~ NA
    )) |>
    unique()


LR_risk_plot_data3 <- LR_risk_plot_data |>
    left_join(LR_risk_plot_data2) |>
    left_join(DE_data |> 
                  select(cell_type = cluster, gene_name, vlmf_t, vlmf_adj.P.Val, vlmf_logFC)) |>
    replace_na(list(mean_mean = 0)) |>
    mutate(interaction = fct_reorder(interaction, mean_mean),
           signif = case_when(vlmf_adj.P.Val < 0.005 ~ "***",
                              vlmf_adj.P.Val < 0.01 ~"**",
                              vlmf_adj.P.Val < 0.05 ~"*",
                              TRUE~"")
    )

LR_risk_plot_data3 |> count(Risk, interaction)

max_abs = max(abs(LR_risk_plot_data3$vlmf_logFC))


LR_risk_plot <- LR_risk_plot_data3 |>
    ggplot(aes(x = cell_type, y = interaction, fill = vlmf_logFC)) +
    facet_grid(Risk~class, scales = "free", space = "free") +
    geom_tile(color = "black") +
    geom_text(aes(label = signif)) +
    scale_fill_gradient2(high = APOE_carrier_colors[["E4+"]], low = APOE_carrier_colors[["E2+"]])  +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom") 

ggsave(LR_risk_plot, filename = here(plot_dir, sprintf("LR_risk_%s.png", opt$datatype)), width = 10, height = 6)


LR_risk_plot_ct <- LR_risk_plot_data3 |>
    ggplot(aes(x = cell_type, y = interaction, fill = cell_type)) +
    facet_grid(Risk~class, scales = "free", space = "free") +
    geom_tile(color = "black") +
    geom_text(aes(label = signif)) +
    scale_fill_manual(values = cluster_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "None") 

ggsave(LR_risk_plot_ct, filename = here(plot_dir, sprintf("LR_risk_ct_%s.png", opt$datatype)), width = 10, height = 6)

#### Bivaraite DEG data ####

LR_DEGS_SpD_bar <- bivariate_SpD_data |>
    mutate(interaction = gsub("\\^", " -> ", pair_id)) |>
    filter(interaction %in% LR_DEGS_plot_data3$interaction) |>
    ggplot(aes(y = interaction, x = mean_score, fill = SpD)) +
    geom_col(position = "dodge") +
    theme_bw()

ggsave(LR_DEGS_SpD_bar, filename = here(plot_dir, "LR_DEGS_bivarite_SpD_bar.png"))

LR_DEGS_SpD_tile<- bivariate_SpD_data |>
    mutate(interaction = gsub("\\^", " -> ", pair_id)) |>
    filter(interaction %in% LR_DEGS_plot_data3$interaction) |>
    ggplot(aes(y = interaction, x = SpD, fill = mean_score)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(LR_DEGS_SpD_tile, filename = here(plot_dir, "LR_DEGS_bivarite_SpD_tile.png"))

# slurmjobs::job_loop(loops = list(datatype = c("sn_broad","sn_fine","Visium"),
#                     create_shell = TRUE,
#                     name = "10_GO_analysis_contrast",
#                     create_script = FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


