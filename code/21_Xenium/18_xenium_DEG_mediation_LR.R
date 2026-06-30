## Louise Huuki-Myers, June 2026
## Explore Xenium DEG results for LR and mediation genes

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("tidyverse")
library("spatialLIBD")
# library("DeconvoBuddies")
# library("scDotPlot")

data_dir <- here("processed-data", "21_Xenium", "18_xenium_DEG_mediation_LR")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "18_xenium_DEG_mediation_LR")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### DEGs - Astro Mediation Genes ####

mediator_outcome_eval <- read_csv(here("processed-data", "21_Xenium", "06_xenium_ALL_probe_eval", "ECR_mediator_outcome_xenium_eval.csv"))

mediator_outcome_eval |> dplyr::count(med_cl)
# med_cl      n
# <chr>   <int>
# 1 Astro.1     5
# 2 Astro.2    35
# 3 Astro.3    18

## cell type level Xenium DEGs

xenium_DEG <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "Xenium_cell_type_anno", "DGE_results_carrier_Xenium_cell_type_anno_wSN.Rds"))

mediator_DEG <- xenium_DEG |>
    inner_join(mediator_outcome_eval |> 
                   select(cluster = med_cl, gene_name = mediator, pair)) |>
    select(cluster, gene_name, vlmf_xenium_logFC, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_xenium_adj.P.Val, validate, pair) |>
    unique()

mediator_DEG |> filter(validate) |> select(-pair) |> unique()

# cluster gene_name vlmf_xenium_logFC vlmf_xenium_t vlmf_xenium_P.Value vlmf_xenium_adj.P.Val validate
# <chr>   <chr>                 <dbl>         <dbl>               <dbl>                 <dbl> <lgl>   
# 1 Astro.1 NPTXR                 0.371          2.71              0.0149                 0.259 TRUE    
# 2 Astro.2 FZD8                  0.337          2.33              0.0330                 0.408 TRUE 


outcome_DEG <- xenium_DEG |>
    inner_join(mediator_outcome_eval |> 
                   select(gene_name = outcome, pair) |>
                   mutate(cluster = "Oligo.3") ) |>
    select(cluster, gene_name, vlmf_xenium_logFC, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_xenium_adj.P.Val, validate, pair)

outcome_DEG |> filter(validate) |> select(-pair) |> unique()

# cluster gene_name vlmf_xenium_logFC vlmf_xenium_t vlmf_xenium_P.Value vlmf_xenium_adj.P.Val validate
# <chr>   <chr>                 <dbl>         <dbl>               <dbl>                 <dbl> <lgl>   
# 1 Oligo.3 SLC17A7               0.221          1.96              0.0669                 0.672 TRUE    
# 2 Oligo.3 TESPA1                0.349          2.00              0.0618                 0.672 TRUE    
# 3 Oligo.3 NPTXR                 0.364          1.92              0.0721                 0.672 TRUE    
# 4 Oligo.3 ENC1                  0.211          1.80              0.0900                 0.726 TRUE  




mediator_outcome_DEG <- mediator_DEG |>
    select(-gene_name) |>
    dplyr::rename(validate_mediator = validate) |>
    rename_with(~str_replace(.x, "xenium", "mediator")) |>
    separate(pair, into = c("mediator", "outcome"), sep = "\\|", remove = FALSE) |> 
    dplyr::rename(cluster_validate = cluster) |>
    left_join(outcome_DEG |>
                  dplyr::rename(outcome = gene_name, 
                                cluster_outcome = cluster,
                                validate_outcome = validate) |>
                  rename_with(~str_replace(.x, "xenium", "outcome"))
              )

mediator_outcome_DEG |> dplyr::count(validate_mediator, validate_outcome)
# validate_mediator validate_outcome     n
# <lgl>             <lgl>            <int>
# 1 FALSE             FALSE               39
# 2 FALSE             TRUE                 4
# 3 TRUE              FALSE               11
# 4 TRUE              TRUE                 4

mediator_outcome_DEG |> dplyr::count(cluster_validate, validate_mediator, validate_outcome)

# cluster_validate validate_mediator validate_outcome     n
# <chr>            <lgl>             <lgl>            <int>
# 1 Astro.1          FALSE             TRUE                 1
# 2 Astro.1          TRUE              FALSE                1
# 3 Astro.1          TRUE              TRUE                 3 *
# 4 Astro.2          FALSE             FALSE               21
# 5 Astro.2          FALSE             TRUE                 3
# 6 Astro.2          TRUE              FALSE               10
# 7 Astro.2          TRUE              TRUE                 1 *
# 8 Astro.3          FALSE             FALSE               18

mediator_outcome_DEG |> 
    filter(validate_mediator | validate_outcome) |> 
    select(pair, cluster_validate, validate_mediator, validate_outcome) |>
    arrange(validate_mediator, validate_outcome)
    

write.csv(mediator_outcome_DEG, file = here(data_dir, "Xenium_mediator_outcome_DEG_results.csv"), row.names = FALSE)
    
 
comapre_stats_scatter <- function(dge_tb, ctb, stat = "t", mX = "vlmf_sn", mY = "vlmf_xenium", FDR_cut_mX = 0.05, pval_cut_mY = 0.1, model_name){
    
    ## define vars
    statX <- paste0(mX, "_", stat)
    statY <- paste0(mY, "_", stat)
    fdrX <- paste0(mX, "_adj.P.Val")
    pvalY <- paste0(mY, "_P.Value")
    # fdrY <- paste0(mY, "_adj.P.Val")
    
    # define colors
    signif_colors <- c("purple", "blue", "red")
    # names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "FDR<", FDR_cut_mY))
    names(signif_colors) <- c("sig_both", paste(mX, "FDR<", FDR_cut_mX) , paste(mY, "P value<", pval_cut_mY))
    
    # make scatter plot
    stat_scatter <- dge_tb |>
        filter(cell_type_broad == ctb) |>
        mutate(DE_class = case_when(!!sym(fdrX) < FDR_cut_mX & !!sym(pvalY) < pval_cut_mY~ "sig_both",
                                    !!sym(fdrX) < FDR_cut_mX ~ paste(mX, "FDR<", FDR_cut_mX),
                                    # !!sym(fdrY) < FDR_cut_mY ~ paste(mY, "FDR<", FDR_cut_mY), ## use pval for xenium data
                                    !!sym(pvalY) < pval_cut_mY ~ paste(mY, "P value<", pval_cut_mY),
                                    TRUE ~ "None")) |>
        ggplot(aes(x = !!sym(statX), y = !!sym(statY), color = DE_class)) +
        geom_point(alpha = 0.5, size = 0.5) +
        geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) +
        geom_abline(linetype = "dashed") +
        scale_color_manual(values = signif_colors) +
        labs(title = model_name, subtitle = paste(mX, "vs.", mY)) +
        facet_wrap(~cluster) +
        theme_bw()
    
    plot_fn = sprintf("%s_%s_stat_scatter_%s_%s-v-%s_%s.png", "Xenium_mediation", stat, model_name, mX, mY, ctb)
    ggsave(stat_scatter, filename = here(plot_dir, plot_fn), height = 7, width = 7)
    
    # return(t_stat_scatter)
}


xenium_DEG_med <- xenium_DEG |>
    inner_join(mediator_outcome |> select(cluster = med_cl, gene_name = mediator)|>
                   unique()) |>
    mutate(cluster = droplevels(cluster))

comapre_stats_scatter(xenium_DEG_out, ctb = "Astro", model_name = "mediator")


xenium_DEG_out <- xenium_DEG |> 
    filter(cluster == "Oligo.3", gene_name %in% mediator_outcome$outcome) |>
    mutate(role = "outcome")

comapre_stats_scatter(xenium_DEG_out, ctb = "Oligo", model_name = "outcome")

xenium_DEG_out |> 
    filter(gene_name == "NPTXR") |>
    select(gene_name, vlmf_t, vlmf_P.Value, vlmf_sn_adj.P.Val, vlmf_sn_logFC, vlmf_sn_t)

head(mediator_outcome)

