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


## add mediator stats
mediator_DE_stats <- mediation_summary |>
    distinct(mediator_datatype, med_cl_test, mediator) |>
    pmap_dfr(function(mediator_datatype, med_cl_test, mediator) {
        DE_data <- load_DE_cached(get_DE_fn(mediator_datatype))
        DE_data |>
            filter(cluster == med_cl_test, gene_name == mediator) |>
            transmute(
                mediator_datatype, med_cl_test, mediator,
                mediatorDE_logFC = vlmf_logFC,
                mediatorDE_P.Value = vlmf_P.Value,
                mediatorDE_adj.P.Val = vlmf_adj.P.Val,
                mediatorDE_t = vlmf_t
            )
    })
mediator_DE_stats |> filter(mediatorDE_P.Value < 0.1)
mediator_DE_stats |> filter(mediator == "SV2B")

mediation_summary <- mediation_summary |> left_join(mediator_DE_stats)

message(Sys.time(), sprintf(" - %d total mediator x outcome rows loaded", nrow(mediation_summary)))

mediation_summary |> filter(base_valid, mediated)

mediated_hits <- mediation_summary |> filter(base_valid, !med_sig)
message(Sys.time(), sprintf(" - %d mediated gene pairs to plot", nrow(mediated_hits)))

mediated_hits |> dplyr::count(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator)

mediation_summary |> 
    filter(base_valid, mediated) |> 
    select( med_cl_test, mediator, outcome_cl, outcome, base_valid, mediated)

mediation_summary |> 
    filter(base_valid, mediated) |> 
    select( med_cl_test, mediator, outcome_cl, outcome, base_valid, mediated, mediatorDE_P.Value) |>
    mutate(mediatorDE = ifelse(mediatorDE_P.Value < 0.1, "*",""))

mediation_summary |> 
    filter(base_valid, !med_sig,  !mediated) |> 
    select( med_cl_test, mediator, outcome_cl, outcome, base_valid, med_sig, med_vec_sig, mediatorDE_P.Value) |>
    arrange(mediator, outcome) |>
    mutate(mediatorDE = ifelse(mediatorDE_P.Value < 0.1, "*",""))

mediated_hits |> select(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator, outcome, Xen_valid)

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
    # med_cl = "Astro.2"
    # med_cl_test = "Astro.2"
    # mediator_datatype = "Xenium_cell_type_anno"
    # outcome_datatype = "Xenium_cell_type_anno"
    # 
    # # med_cl_test = "Astro.2_APOE_high"
    # # mediator_datatype = "Xenium_Oligo.3_Astro"
    # # outcome_datatype = "Xenium_Oligo.3_Astro"
    # 
    # mediator = "SV2B"
    # outcome = "ENC1"
    # outcome_cl = "APOE_high_nnA_Astro.2"
    # outcome_cl = "Oligo.3"
    # 
    message(Sys.time(), sprintf(" - plotting mediator=%s (%s) -> %s (%s)", mediator, med_cl_test, outcome, outcome_cl))
    # 
    sce_out <- load_pb_cached(get_pb_fn(outcome_datatype))
    sce_med <- load_pb_cached(get_pb_fn(mediator_datatype))
    mediator_stats <- load_DE_cached(get_DE_fn(mediator_datatype))

    med_cl_test %in% mediator_stats$cluster
    
    outcome_cl %in% sce_out[["registration_variable"]]
    
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
        )), height = 4, width = 8)
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
