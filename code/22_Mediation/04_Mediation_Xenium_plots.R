## Louise Huuki-Myers, July 2026
## Adapted from 04_DEG_boxplots.R to plot the mediated gene pairs found by
## 03_Mediation_Xenium_pairs.R, using plot_DEG_mediated_express() (the
## two-panel unadjusted-vs-adjusted-for-mediator version of
## plot_DEG_express()).

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("DeconvoBuddies")
library("SingleCellExperiment")
library("patchwork")

source(here("code", "utils", "plot_DEG_express.R"))
source(here("code", "22_Mediation", "plot_DEG_mediated_express.R"))
source(here("code", "22_Mediation", "xenium_datatype_lookup.R"))

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

plot_dir <- here("plots", "22_Mediation", "04_Mediation_Xenium_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load xenium mediation data ####
mediation_summary_fn <- here("processed-data", "22_Mediation", "03_Mediation_Xenium", "mediation_summary-all_scenarios_Pval0.10.csv")
mediation_summary <- read_csv(mediation_summary_fn)

mediation_summary |> dplyr::count(mediator, outcome) |> arrange(-n)
# A tibble: 58 × 3
# mediator outcome     n
# <chr>    <chr>   <int>
# 1 ABCA8    MBP        10
# 2 FZD8     CABLES1    10
# 3 FZD8     CPNE4      10
# 4 FZD8     GAD1       10
# 5 FZD8     GPM6A      10
# ...

mediation_summary |> 
    filter(mediator == "NPTXR", outcome == "GAD1") |> 
    select(mediator, outcome, outcome_cl, P.Value_base, P.Value_med)

## add mediator stats
mediator_DE_stats <- mediation_summary |>
    distinct(mediator_datatype, med_cl_test, mediator) |>
    pmap_dfr(function(mediator_datatype, med_cl_test, mediator) {
        DE_data <- load_DE_cached(get_DE_fn(mediator_datatype))
        # SpX DE tables store the full spatial cluster label in cluster_xSpD;
        # cluster is stripped to the base cell type (e.g. "Astro.2" vs
        # "Astro.2_L6"). med_cl_test uses the full label, so remap for SpX.
        if (grepl("_SpX$", mediator_datatype)) {
            DE_data <- DE_data |> mutate(cluster = cluster_xSpD)
        }
        DE_data |>
            filter(cluster == med_cl_test, gene_name == mediator) |>
            transmute(
                mediator_datatype, 
                med_cl_test, 
                mediator,
                mediatorDE_logFC = vlmf_logFC,
                mediatorDE_P.Value = vlmf_P.Value,
                mediatorDE_adj.P.Val = vlmf_adj.P.Val,
                mediatorDE_t = vlmf_t
            )
    })

mediator_DE_stats |> dplyr::count(is.na(mediatorDE_P.Value))
mediator_DE_stats |> dplyr::count(mediator_datatype, med_cl_test)

# 4 validated mediators
mediator_DE_stats |> filter(mediatorDE_P.Value < 0.1) |> dplyr::count(mediator)
# mediator     n
# <chr>    <int>
# 1 FZD8         5
# 2 NPTXR        4
# 3 ST18         1
# 4 SV2B         1

mediator_DE_stats |> filter(mediator == "SV2B")

mediation_summary <- mediation_summary |> left_join(mediator_DE_stats)

message(Sys.time(), sprintf(" - %d total mediator x outcome rows loaded", nrow(mediation_summary)))
# 408 total mediator x outcome rows loaded

mediation_summary |> filter(base_valid, mediated)

mediated_hits <- mediation_summary |> filter(base_valid, !med_sig)

mediated_hits |> dplyr::count(mediator, outcome)

# mediator outcome     n
# <chr>    <chr>   <int>
# 1 FZD8     CABLES1     2
# 2 FZD8     CPNE4       1
# 3 FZD8     GAD1        1
# 4 FZD8     NTRK3       1
# 5 FZD8     TESPA1      4
# 6 NPTXR    ENC1        1
# 7 NPTXR    GAD1        1
# 8 NPTXR    NPTXR       2
# 9 NPTXR    SLC17A7     1
# 10 SV2B    ERBB3       1

message(Sys.time(), sprintf(" - %d mediated gene pairs to plot", nrow(mediated_hits)))
# 15 mediated gene pairs to plot
mediated_hits |> dplyr::count(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator, outcome) 

mediation_summary |> 
    filter(base_valid, mediated) |> 
    arrange(mediator) |>
    select( med_cl_test, mediator,  P.Value_med, outcome_cl, outcome, base_valid, mediated)

# med_cl_test       mediator outcome_cl            outcome base_valid mediated
# <chr>             <chr>    <chr>                 <chr>   <lgl>      <lgl>   
# 1 Astro.2_xL5       FZD8     Oligo.3_xL5           CABLES1 TRUE       TRUE    
# 2 Astro.1           NPTXR    Oligo.3               NPTXR   TRUE       TRUE    
# 3 Astro.1_xLD_glia  NPTXR    Oligo.3_xLD_glia      NPTXR   TRUE       TRUE    
# 4 Astro.2_APOE_high SV2B     APOE_high_nnA_Astro.2 ERBB3   TRUE       TRUE 

mediation_summary |> 
    filter(base_valid, mediated) |> 
    select( med_cl_test, mediator, outcome_cl, outcome, base_valid, mediated, mediatorDE_P.Value) |>
    mutate(mediatorDE = ifelse(mediatorDE_P.Value < 0.1, "*",""))

mediation_summary |> 
    filter(base_valid, !med_sig,  !mediated) |> 
    select( med_cl_test, mediator, outcome_cl, outcome, base_valid, med_sig, med_vec_sig, mediatorDE_P.Value) |>
    arrange(mediator, outcome) |>
    mutate(mediatorDE = ifelse(mediatorDE_P.Value < 0.1, "*",""))

mediated_hits |> select(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator, outcome, base_valid)

write_csv(mediated_hits, file = here("processed-data", "22_Mediation", "03_Mediation_Xenium", "Xenium_mediation_hits.csv"))

#### Test ####

sce_pb_test <- readRDS(here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno.RDS"))
sce_pb_test$APOE_carrier_syn <- gsub("\\+", "", sce_pb_test$APOE_carrier)

## mediator stats 
mediated_hits_test <- mediated_hits |> filter(med_cl_test == "Astro.1")

med_plot_test <- plot_DEG_mediated_express(
    sce = sce_pb_test,
    sce_mediator = sce_pb_test,
    stats = mediation_summary,
    clus = "Oligo.3",
    med_clus = "Astro.1",
    mediator_gene = "NPTXR",
    gene = "GAD1",
    gene_col = "gene_name",
    cluster_col = "registration_variable",
    med_cluster_col = "registration_variable",
    category_col = "APOE_carrier",
    mod = ~APOE_carrier_syn + Age + Anc_Afr,
    cleanY_P = 2,
    color_pal = APOE_carrier_colors,
    signif_stat = "P.Value",
    signif_thr = 0.10
)

ggsave(med_plot_test, filename = here(plot_dir, "med_plot_test.png"))

# DE_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds"))
DE_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "Xenium_cell_type_anno", "DGE_results_carrier_Xenium_cell_type_anno.Rds"))

med_plot_test2 <- plot_DEG_mediated_express(
    sce = sce_pb_test,
    sce_mediator = sce_pb_test,
    stats = mediation_summary,
    clus = "Oligo.3",
    med_clus = "Astro.1",
    mediator_gene = "NPTXR",
    gene = "GAD1",
    mod = ~APOE_carrier_syn + Age + Anc_Afr,
    cleanY_P = 2,
    color_pal = APOE_carrier_colors,
    plot_points = TRUE,
    plot_mediator_panel = TRUE,
    mediator_stats = DE_data,
    mediator_stats_cluster_col = "cluster",
    mediator_pval_col = "vlmf_P.Value",
    mediator_fc_col = "vlmf_logFC",
    mediator_mod = ~APOE_carrier_syn + Age + Anc_Afr,
    mediator_cleanY_P = 2,
    signif_stat = "P.Value",
    signif_thr = 0.10
)

ggsave(med_plot_test2, filename = here(plot_dir, "med_plot_test2.png"))



## pb_lookup / get_pb_fn / load_pb_cached and DE_lookup / get_DE_fn /
## load_DE_cached now live in xenium_datatype_lookup.R (sourced above),
## shared with 03_Mediation_Xenium_pairs.R so the two scripts can't drift
## apart on which file a given datatype key points to.

#### model for cleaning outcome expression -- matches 06_Clusterwise_voomLmFit_Xenium.R ####
mediation_mod <- ~APOE_carrier_syn + Age + Anc_Afr
mediation_cleanY_P <- 2  ## keeps (Intercept) + APOE_carrier_E4, regresses out Age + Anc_Afr


#### Loop over each (med_cl, mediator) group and plot ####

mediated_hits_select <- mediated_hits |> select(med_cl, med_cl_test, mediator_datatype, outcome_datatype, mediator, outcome, outcome_cl)

pdf(here(plot_dir, "mediation_boxplots_Xenium_pairs.pdf"), width = 10, height = 4)

pwalk(mediated_hits_select, function(med_cl, med_cl_test, mediator_datatype, outcome_datatype, mediator, outcome, outcome_cl) {
    
    
    # ## test
    # 14 Astro.2 Astro.2_L6        Xenium_cell_type_anno_SpX Xenium_cell_type_anno_SpX FZD8     CABLES1 Oligo.3_L6
    # med_cl = "Astro.2"
    # med_cl_test = "Astro.2_L6"
    # mediator_datatype = "Xenium_cell_type_anno_SpX"
    # outcome_datatype = "Xenium_cell_type_anno_SpX"
    # 
    # # med_cl_test = "Astro.2_APOE_high"
    # # mediator_datatype = "Xenium_Oligo.3_Astro"
    # # outcome_datatype = "Xenium_Oligo.3_Astro"
    # 
    # mediator = "FZD8"
    # outcome = "CABLES1"
    # outcome_cl = "Oligo.3_L6"
    # outcome_cl = "Oligo.3"
    # 
    message(Sys.time(), sprintf(" - plotting mediator=%s (%s) -> %s (%s)", mediator, med_cl_test, outcome, outcome_cl))
    # 
    sce_out <- load_pb_cached(get_pb_fn(outcome_datatype))
    sce_med <- load_pb_cached(get_pb_fn(mediator_datatype))
    mediator_stats <- load_DE_cached(get_DE_fn(mediator_datatype))
    
    if(grepl("_SpX$", mediator_datatype)) mediator_stats <- mediator_stats |> mutate(cluster = cluster_xSpD)

    # med_cl_test %in% mediator_stats$cluster
    # outcome_cl %in% sce_out[["registration_variable"]]
    
    p <- tryCatch(
        plot_DEG_mediated_express(
            sce = sce_out,
            sce_mediator = sce_med,
            stats = mediated_hits,
            clus = outcome_cl,
            med_clus = med_cl_test,
            mediator_gene = mediator,
            gene = outcome,
            gene_col = "gene_name",
            cluster_col = "registration_variable",
            med_cluster_col = "registration_variable",
            mod = mediation_mod,
            cleanY_P = mediation_cleanY_P,
            color_pal = APOE_carrier_colors,
            plot_points = TRUE,
            plot_mediator_panel = TRUE,
            mediator_stats = mediator_stats,
            mediator_stats_cluster_col = "cluster",
            signif_stat = "P.Value",
            signif_thr = 0.10
        ),
        error = function(e) {
            warning(sprintf("Failed to plot mediator=%s med_cl_test=%s: %s", mediator, med_cl_test, conditionMessage(e)))
            NULL
        }
    )
    
    if (!is.null(p)) {
        print(p + patchwork::plot_annotation(title = sprintf("%s (%s) -> %s", mediator, med_cl, outcome_cl)))
        ggsave(p, filename = here(plot_dir, sprintf(
            "mediation_Xen_%s_%s_%s_%s.png", gsub("\\.", "", med_cl), mediator, gsub("\\.", "", outcome_cl), outcome
        )), height = 3.5, width = 7)
    }
})

dev.off()



#### Loop over each (scenario, mediator) group and plot all its mediated outcome genes together ####

plot_groups <- mediated_hits |>
    distinct(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator)

message(Sys.time(), sprintf(" - %d scenario x mediator groups to plot", nrow(plot_groups)))

pdf(here(plot_dir, "mediation_boxplots_all_scenarios.pdf"), width = 8, height = 4)

pwalk(plot_groups, function(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator) {

    message(Sys.time(), sprintf(" - plotting mediator=%s (%s) -> outcome_cl=%s (%s)",
                                 mediator, med_cl_test, outcome_cl, outcome_datatype))

    sce_out <- load_pb_cached(get_pb_fn(outcome_datatype))
    sce_med <- load_pb_cached(get_pb_fn(mediator_datatype))

    genes_this_group <- mediated_hits |>
        filter(mediator_datatype == !!mediator_datatype, 
               outcome_datatype == !!outcome_datatype,
               outcome_cl == !!outcome_cl, 
               med_cl_test == !!med_cl_test, 
               mediator == !!mediator) |>
        pull(outcome)

    p <- tryCatch(
        plot_DEG_mediated_express(
            sce = sce_out,
            sce_mediator = sce_med,
            stats = mediated_hits,
            clus = outcome_cl,
            med_clus = med_cl_test,
            mediator_gene = mediator,
            gene = genes_this_group,
            gene_col = "gene_name",
            cluster_col = "registration_variable",
            med_cluster_col = "registration_variable",
            category_col = "APOE_carrier",
            mod = mediation_mod,
            cleanY_P = mediation_cleanY_P,
            color_pal = APOE_carrier_colors,
            plot_points = TRUE
        ),
        error = function(e) {
            warning(sprintf("Failed to plot mediator=%s outcome_cl=%s: %s", mediator, outcome_cl, conditionMessage(e)))
            NULL
        }
    )

    if (!is.null(p)) {
        print(p + patchwork::plot_annotation(
            title = sprintf("%s | mediator: %s (%s)", outcome_cl, mediator, med_cl_test),
            subtitle = sprintf("datatype: %s (outcome) / %s (mediator)", outcome_datatype, mediator_datatype)
        ))

        ## also save an individual PNG per group, matching the
        ## single-gene ggsave pattern from 04_DEG_boxplots.R
        ggsave(p, filename = here(plot_dir, sprintf(
            "mediation_%s_%s_%s_%s.png",
            gsub("\\.", "", outcome_cl), gsub("\\.", "", med_cl_test), mediator, outcome_datatype
        )), height = 4, width = 3 + 2 * length(genes_this_group))
    }
})

dev.off()

message(Sys.time(), sprintf(" - Done: plotted %d scenario x mediator groups covering %d mediated gene pairs",
                             nrow(plot_groups), nrow(mediated_hits)))

# slurmjobs::job_single('04_Mediation_Xenium_plots', create_shell = TRUE, memory = '10G',
#                       command = "Rscript 04_Mediation_Xenium_plots.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
