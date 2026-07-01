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

#### Liana LR Data ####

LR_pair_df = read_csv(here('processed-data', '24_xenium_liana', '04_global_score_heatmap', 'unique_ligand_receptor_pairs.csv'), show_col_types = FALSE) |>
    mutate(LR_pair = paste0(ligand_complex, "|", receptor_complex))

LR_pair_df |> filter("APOE" == ligand_complex)
LR_pair_df |> filter("APP" == receptor_complex)

LR_pair_df |> filter(LR_pair %in% mediator_outcome_eval$pair) ## no overlapping pairs

LR_pair_df |> filter(ligand_complex %in% mediator_outcome_eval$mediator) # no overlap

(mediator_LR_overlap <- LR_pair_df |> 
    filter(receptor_complex %in% mediator_outcome_eval$outcome) |>
    left_join(mediator_outcome_eval |> 
                  select(receptor_complex = outcome, mediator_pair = pair, mediator),
              relationship = "many-to-many"))

# ligand_complex receptor_complex LR_pair      mediator_pair mediator
# <chr>          <chr>            <chr>        <chr>         <chr>   
# 1 NXPH2          NRXN1            NXPH2|NRXN1  FZD8|NRXN1    FZD8    
# 2 NXPH2          NRXN1            NXPH2|NRXN1  NETO1|NRXN1   NETO1   
# 3 NXPH2          NRXN1            NXPH2|NRXN1  SV2B|NRXN1    SV2B    
# 4 NXPH2          NRXN1            NXPH2|NRXN1  CHRM3|NRXN1   CHRM3   
# 5 NLGN1          NRXN1            NLGN1|NRXN1  FZD8|NRXN1    FZD8    
# 6 NLGN1          NRXN1            NLGN1|NRXN1  NETO1|NRXN1   NETO1   
# 7 NLGN1          NRXN1            NLGN1|NRXN1  SV2B|NRXN1    SV2B    
# 8 NLGN1          NRXN1            NLGN1|NRXN1  CHRM3|NRXN1   CHRM3   

# 9 NRG3           ERBB3            NRG3|ERBB3   SV2B|ERBB3    SV2B    
# 10 S100A4         ERBB3            S100A4|ERBB3 SV2B|ERBB3    SV2B 


xenium_DEG <- xenium_DEG |>
    mutate(ligand = gene_name %in% LR_pair_df$ligand_complex,
           receptor = gene_name %in% LR_pair_df$receptor_complex,
           mediator = gene_name %in% mediator_outcome_eval$mediator,
           outcome = gene_name %in% mediator_outcome_eval$outcome)

xenium_DEG |> dplyr::count(ligand, receptor)

xenium_DEG |> 
    filter(cluster == "Oligo.3", receptor) |>
    arrange(vlmf_xenium_P.Value) |>
    filter(signif_sn | signif_xenium) |>
    select(gene_name, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_sn_adj.P.Val, vlmf_sn_logFC, vlmf_sn_t, signif_sn, signif_xenium, validate, outcome)

# gene_name vlmf_xenium_t vlmf_xenium_P.Value vlmf_sn_adj.P.Val vlmf_sn_logFC vlmf_sn_t signif_sn signif_xenium validate
# <chr>             <dbl>               <dbl>             <dbl>         <dbl>     <dbl> <lgl>     <lgl>         <lgl>   
# 1 NLGN1           -2.54                0.0210            0.0739        -0.538     -2.74 FALSE     TRUE          FALSE   
# 2 CAV1            -2.29                0.0350            0.385         -0.572     -1.38 FALSE     TRUE          FALSE   
# 3 PTPRD           -2.22                0.0401            0.0245        -0.493     -3.77 TRUE      TRUE          TRUE    **
# 4 ERBB3           -0.749               0.464             0.0429        -0.601     -3.20 TRUE      FALSE         FALSE   
# 5 NRXN3           -0.521               0.609             0.0191        -0.684     -4.06 TRUE      FALSE         FALSE   
# 6 NRXN1           -0.422               0.678             0.0398         1.15       3.27 TRUE      FALSE         FALSE   
# 7 TLR2             0.0164              0.987             0.0394         1.56       3.28 TRUE      FALSE         FALSE 

xenium_DEG |> 
    filter(cluster == "Oligo.3", receptor,
           gene_name %in% c("APP", "TREM2", "FZD8")) |>
    select(gene_name, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_sn_adj.P.Val, vlmf_sn_logFC, vlmf_sn_t, signif_sn, signif_xenium, validate, outcome)

# gene_name vlmf_xenium_t vlmf_xenium_P.Value vlmf_sn_adj.P.Val vlmf_sn_logFC vlmf_sn_t signif_sn signif_xenium validate outcome
# <chr>             <dbl>               <dbl>             <dbl>         <dbl>     <dbl> <lgl>     <lgl>         <lgl>    <lgl>  
# 1 APP               1.18                0.254             0.106        -0.289     -2.45 FALSE     FALSE         FALSE    FALSE  
# 2 TREM2            -0.116               0.909            NA            NA         NA    NA        FALSE         FALSE    FALSE 

xenium_DEG |> 
    filter(grepl("Astro", cluster), ligand) |>
    arrange(gene_name, vlmf_xenium_P.Value) |>
    filter(signif_sn | signif_xenium) |>
    select(cluster, gene_name, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_sn_adj.P.Val, vlmf_sn_logFC, vlmf_sn_t, signif_sn, signif_xenium, validate, mediator)

xenium_DEG |> 
    filter(grepl("Astro", cluster),
           gene_name %in% c("APOE", "SPON1")) |>
    arrange(gene_name, vlmf_xenium_P.Value) |>
    select(cluster, gene_name, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_sn_adj.P.Val, vlmf_sn_logFC, vlmf_sn_t, signif_sn, signif_xenium, validate, mediator)

# cluster gene_name vlmf_xenium_t vlmf_xenium_P.Value vlmf_sn_adj.P.Val vlmf_sn_logFC vlmf_sn_t signif_sn signif_xenium validate
# <fct>   <chr>             <dbl>               <dbl>             <dbl>         <dbl>     <dbl> <lgl>     <lgl>         <lgl>   
# 1 Astro.2 APOE              0.898              0.382              0.630      0.330      1.43    FALSE     FALSE         FALSE   
# 2 Astro.1 APOE              0.550              0.589              0.970     -0.0412    -0.202   FALSE     FALSE         FALSE   
# 3 Astro.5 APOE              0.549              0.590              0.646      0.411      1.50    FALSE     FALSE         FALSE   
# 4 Astro.3 APOE              0.461              0.650              0.834      0.195      0.917   FALSE     FALSE         FALSE   
# 5 Astro.4 APOE             -0.104              0.918              0.664     -0.491     -1.97    FALSE     FALSE         FALSE 

# 6 Astro.3 SPON1            -2.76               0.0134             0.924     -0.0660    -0.499   FALSE     TRUE          FALSE   
# 7 Astro.4 SPON1            -2.62               0.0177             0.646     -0.615     -2.04    FALSE     TRUE          FALSE   # SPON1 nominaly downregulated in Astrocytes
# 8 Astro.1 SPON1            -2.01               0.0606             0.999      0.000320   0.00310 FALSE     TRUE          FALSE   
# 9 Astro.5 SPON1            -1.36               0.191              0.958      0.0420     0.252   FALSE     FALSE         FALSE   
# 10 Astro.2 SPON1            -1.16               0.263              0.836     -0.0996    -0.768   FALSE     FALSE         FALSE


xenium_DEG |> 
    filter(
        grepl("Astro", cluster),
           gene_name %in% mediator_LR_overlap$mediator,
           validate) |>
    arrange(gene_name, vlmf_xenium_P.Value) |>
    select(cluster, gene_name, vlmf_xenium_t, vlmf_xenium_P.Value, vlmf_sn_adj.P.Val, vlmf_sn_logFC, vlmf_sn_t, signif_sn, signif_xenium, validate, mediator)

# cluster gene_name vlmf_xenium_t vlmf_xenium_P.Value vlmf_sn_adj.P.Val vlmf_sn_logFC vlmf_sn_t signif_sn signif_xenium validate mediator
# <fct>   <chr>             <dbl>               <dbl>             <dbl>         <dbl>     <dbl> <lgl>     <lgl>         <lgl>    <lgl>   
#1 Astro.2 FZD8               2.33              0.0330            0.0261          1.01      4.63 TRUE      TRUE          TRUE     TRUE   

mediator_LR_overlap |> filter(mediator == "CHRM3")
mediator_LR_overlap |> filter(mediator == "FZD8")

mediator_outcome_eval |> filter(mediator == "FZD8")

#### Check mediation genes in LR ####


mediator_outcome_DEG |>
    filter(mediator %in% pr)


