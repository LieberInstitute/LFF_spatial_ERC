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

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

plot_dir <- here("plots", "22_Mediation", "04_Mediation_Xenium_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load mediation data ####
mediation_summary_fn <- here("processed-data", "22_Mediation", "03_Mediation_Xenium", "mediation_summary-all_scenarios_Pval0.10.csv")
mediation_summary <- read_csv(mediation_summary_fn)

message(Sys.time(), sprintf(" - %d total mediator x outcome rows loaded", nrow(mediation_summary)))

mediated_hits <- mediation_summary |> filter(Xen_mediated)
message(Sys.time(), sprintf(" - %d mediated gene pairs to plot", nrow(mediated_hits)))

if (nrow(mediated_hits) < 1) stop("No mediated=TRUE rows found in mediation_summary -- nothing to plot.")

mediated_hits |> dplyr::count(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator)

mediated_hits |> select(mediator_datatype, outcome_datatype, outcome_cl, med_cl_test, mediator, outcome, Xen_valid)

#### Test ####

sce_pb_test <- readRDS(here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno.RDS"))
sce_pb_test$APOE_carrier_syn <- gsub("\\+", "", sce_pb_test$APOE_carrier)

Pvalthr <- 0.10

## mediator stats 
# de_stats <- readRDS(here("processed-data", "13_compile_DGE",))

med_plot_test <- plot_DEG_mediated_express(
    sce = sce_pb_test,
    sce_mediator = sce_pb_test,
    stats = mediation_summary,
    clus = "Oligo.3",
    med_clus = "Astro.1",
    mediator = "NPTXR",
    gene = "ENC1",
    gene_col = "gene_name",
    cluster_col = "registration_variable",
    med_cluster_col = "registration_variable",
    category_col = "APOE_carrier",
    mod = ~APOE_carrier_syn + Age + Anc_Afr,
    cleanY_P = 2,
    color_pal = APOE_carrier_colors
)

ggsave(med_plot_test, filename = here(plot_dir, "med_plot_test.png"))



## pb_fn resolution + caching -- same lookup as 03_Mediation_Xenium_pairs.R.
## Kept in sync manually for now; consider moving this into a shared
## utils/pb_lookup.R sourced by both scripts so the two can't drift apart.

pb_dir <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep")
pb_lookup <- tribble(
    ~datatype_key,                               ~pb_fn,
    "Xenium_cell_type_anno",                     file.path(pb_dir, "spe_xenium_pseudo_DGE-cell_type_anno.RDS"),
    "Xenium_cell_type_anno_SpX",                 file.path(pb_dir, "spe_xenium_pseudo_DGE-cell_type_anno_SpX.RDS"),
    "Xenium_Oligo.3_Astro",                      file.path(pb_dir, "spe_xenium_pseudo_DGE-Oligo.3_Astro.RDS"),
    "Xenium_Oligo.3_Astro_SpX",                  file.path(pb_dir, "spe_xenium_pseudo_DGE-Oligo.3_Astro_SpX.RDS")
)

get_pb_fn <- function(datatype_key) {
    hit <- pb_lookup$pb_fn[pb_lookup$datatype_key == datatype_key]
    if (length(hit) != 1) stop(sprintf("Unrecognized datatype key '%s' -- add it to pb_lookup.", datatype_key))
    hit
}

.pb_cache <- new.env()
load_pb_cached <- function(fn) {
    if (!exists(fn, envir = .pb_cache, inherits = FALSE)) {
        message(Sys.time(), " - loading pseudobulk: ", basename(fn))
        sce <- readRDS(fn)
        sce$APOE_carrier_syn <- gsub("\\+", "", sce$APOE_carrier)
        assign(fn, sce, envir = .pb_cache)
    }
    get(fn, envir = .pb_cache, inherits = FALSE)
}

#### model for cleaning outcome expression -- matches 06_Clusterwise_voomLmFit_Xenium.R ####
mediation_mod <- ~APOE_carrier_syn + Age + Anc_Afr
mediation_cleanY_P <- 2  ## keeps (Intercept) + APOE_carrier_E4, regresses out Age + Anc_Afr

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
        filter(mediator_datatype == !!mediator_datatype, outcome_datatype == !!outcome_datatype,
               outcome_cl == !!outcome_cl, med_cl_test == !!med_cl_test, mediator == !!mediator) |>
        pull(outcome)

    p <- tryCatch(
        plot_DEG_mediated_express(
            sce = sce_out,
            sce_mediator = sce_med,
            stats = mediated_hits,
            clus = outcome_cl,
            med_clus = med_cl_test,
            mediator = mediator,
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
