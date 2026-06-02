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

xenium_DEG <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "Xenium", "DGE_results_carrier_Xenium_wSN.Rds"))

# list.files(here("processed-data","22_Mediation","out-erc_astro"))

mediation_summary <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediation_vs_baseline_summary.tsv.gz"))

mediator_tab <- mediation_summary |> 
    filter(mediated_n > 0) |> 
    select(mediation_run, mediated_n) |> 
    separate(mediation_run, into = c("cell_type", "ensemblID", "gene_name"), sep = "\\|") |>
    mutate(in_base = gene_name %in% xenium_DEG$gene_name)

head(mediator_tab)

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(mediator_in_DEG = mediator %in% xenium_DEG$gene_name,
           outcome_in_DEG = outcome %in% xenium_DEG$gene_name,
           pair = paste0(mediator, "|", outcome))

mediator_outcome |> filter(mediator_in_DEG, outcome_in_DEG)

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

