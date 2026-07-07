## Louise Huuki-Myers June 2026 (base voomLmFit script)
## Mediation step added to validate the snRNA-seq mediation screen
## (01_mediation_screening.R) in Xenium spatial pseudobulk data.
##
## Restructured from 03_Mediation_Xenium.R (single hardcoded scenario) to
## loop through every row of Oligo3_Astro_mediator_outcome_pairs.csv.
## Each row is a "scenario": a (mediator pseudobulk object, outcome
## pseudobulk object, outcome cluster, mediator cluster label) combo.
## Within each scenario, candidate mediator GENES come from
## ECR_mediator_outcome_xenium_eval.csv, filtered to that scenario's
## med_cl -- same discovery-anchored candidate list logic as before, just
## replayed across every pseudobulk scenario instead of one hardcoded run.

#### Set up ####

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("here")
library("getopt")
library("tidyverse")
library("sessioninfo")
library("qs2")

# mopt <- matrix(
#     c("mediatorPval", "m", "1", "double", "P-value threshold for calling a gene mediated (default: 0.10 for validation)"),
#     ncol = 5, byrow = TRUE
# )
# opt <- getopt(mopt)
# if (is.null(opt$mediatorPval) || is.na(opt$mediatorPval)) opt$mediatorPval <- 0.10

Pvalthr <- 0.10

batch <- "chip"
baseline_formula_str <- "APOE_carrier_syn + Age + Anc_Afr" ## must match 06_Clusterwise_voomLmFit_Xenium.R exactly

#### Define functions ####
## pb_fn resolution -- covers BOTH the short datatype codes used by
## 06_Clusterwise_voomLmFit_Xenium.R's --datatype option AND the raw file
## basenames that show up in the pairs CSV's mediator_datatype column
## (e.g. "spe_xenium_pseudo_DGE-Oligo.3_Astro" is the same file as the
## datatype code "Xenium_Oligo.3_Astro" -- worth standardizing the pairs
## CSV to one convention, but both resolve correctly here either way).

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

## DE_data_dir (where 06_Clusterwise_voomLmFit_Xenium.R saved baseline
## topTables) is always keyed by the CLEAN datatype code, never the raw
## basename variant -- confirm this holds for any new datatype values you
## add to the pairs CSV.
get_DE_data_dir <- function(datatype_key) {
    here("processed-data", "12_voomLmFit", "06_Clusterwise_voomLmFit_Xenium", sprintf("vlmf_%s", datatype_key))
}

out_dir <- here("processed-data", "22_Mediation", "03_Mediation_Xenium")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


## Cached loaders -- pb objects are large and several scenario rows share
## the same underlying file (e.g. rows 1-3 of the pairs CSV all point to
## the same mediator AND outcome pb file), so load each file only once.

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

.dge_cache <- new.env()
get_filtered_dge <- function(sce_pb, clus, formula_str) {
    dge <- sce_pb[, sce_pb$registration_variable == clus]
    des <- model.matrix(as.formula(paste("~", formula_str)), data = colData(dge))
    dge <- edgeR::calcNormFactors(dge)
    keep <- edgeR::filterByExpr.DGEList(dge, design = as.data.frame(des))
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    edgeR::calcNormFactors(dge)
}
get_filtered_dge_cached <- function(pb_fn, clus, formula_str) {
    key <- paste(pb_fn, clus, sep = "||")
    if (!exists(key, envir = .dge_cache, inherits = FALSE)) {
        sce_pb <- load_pb_cached(pb_fn)
        assign(key, get_filtered_dge(sce_pb, clus, formula_str), envir = .dge_cache)
    }
    get(key, envir = .dge_cache, inherits = FALSE)
}

.baseline_cache <- new.env()
get_baseline_tt_cached <- function(outcome_datatype, out_clus) {
    key <- paste(outcome_datatype, out_clus, sep = "||")
    if (!exists(key, envir = .baseline_cache, inherits = FALSE)) {
        fn <- here(get_DE_data_dir(outcome_datatype), sprintf("voomLmFit_%s_%s.rds", outcome_datatype, out_clus))
        if (!file.exists(fn)) {
            stop(sprintf("Baseline output not found: %s -- run 06_Clusterwise_voomLmFit_Xenium.R --datatype %s first.",
                         fn, outcome_datatype))
        }
        assign(key, readRDS(fn), envir = .baseline_cache)
    }
    get(key, envir = .baseline_cache, inherits = FALSE)
}


## Core fitting / classification functions (unchanged logic from
## 03_Mediation_Xenium.R -- gene_name is the shared identifier across
## baseline/mediation topTables and pb object rownames for Xenium's
## targeted panel).

get_mediator_vector <- function(sce_pb, med_cluster, gene_id) {
    stopifnot(gene_id %in% rownames(sce_pb))
    stopifnot(med_cluster %in% sce_pb$registration_variable)

    clus_sce <- sce_pb[gene_id, sce_pb$registration_variable == med_cluster]
    brnum <- as.character(colData(clus_sce)$BrNum)
    stopifnot(!anyDuplicated(brnum))

    mediator_vector <- as.numeric(assay(clus_sce, "logcounts"))
    names(mediator_vector) <- brnum
    mediator_vector
}

fit_mediation_de <- function(dge_sub, baseline_formula_str, batch_col, med_vec) {
    sample_df <- dge_sub$samples
    sample_df$med_vec <- med_vec

    des <- model.matrix(as.formula(paste("~", baseline_formula_str, "+ med_vec")), data = sample_df)
    des <- as.data.frame(des)
    colnames(des) <- gsub(colnames(des), pattern = "_syn", replacement = "_")

    v.swt <- voomLmFit(dge_sub, design = des, block = as.factor(sample_df[[batch_col]]),
                        adaptive.span = TRUE, sample.weights = TRUE)
    v.swt.fit.e <- eBayes(v.swt)

    tt <- topTable(v.swt.fit.e, coef = "APOE_carrier_E4", number = Inf, adjust.method = "BH") |>
        as_tibble(rownames = "outcome_gene_id")
    ttM <- topTable(v.swt.fit.e, coef = "med_vec", number = Inf, adjust.method = "BH") |>
        as_tibble(rownames = "outcome_gene_id")

    list(tt = tt, ttM = ttM)
}

## discovery-anchored validation details for one (med_cluster, mediator gene)
## pair -- restricted to the specific outcome genes that pair was
## significant for in the snRNA-seq screen (mediation_eval), with
## directional concordance against the discovery t-statistic.
mediation_validation_details <- function(mediation_eval, med_cluster, gene_name,
                                          tt_baseline, tt_mediation, ttM_mediation, Pvalthr = 0.10) {
    mediation_eval |>
        filter(med_cl == med_cluster, mediator == gene_name) |>
        select(med_cl:t_med) |>
        rename_with(~paste0(.x, "_SN_outcome"), fdr:t_med) |>
        left_join(tt_baseline   |> select(outcome = gene_name, P.Value_Xen = P.Value,    t_Xen = t),    by = join_by(outcome)) |>
        left_join(tt_mediation  |> select(outcome = gene_name, P.Value_Xen_carrier = P.Value, t_Xen_carrier = t), by = join_by(outcome)) |>
        left_join(ttM_mediation |> select(outcome = gene_name, P.Value_Xen_med = P.Value,     t_Xen_med = t),     by = join_by(outcome)) |>
        mutate(
            Xen_sig    = P.Value_Xen < Pvalthr,
            Xen_dir    = sign(t_SN) == sign(t_Xen),
            Xen_valid  = Xen_sig & Xen_dir,
            Xen_carrier_sig = P.Value_Xen_carrier < Pvalthr,
            Xen_carrier_dir = sign(t_Xen_carrier) == sign(t_Xen),
            Xen_med_sig     = P.Value_Xen_med < Pvalthr,
            Xen_mediated    = !Xen_carrier_sig & Xen_med_sig
        )
}


#### Load the two driver tables ####
##  - pairs_tbl: which pseudobulk scenario to run (outer loop)
##  - mediation_eval: which specific mediator/outcome gene pairs to test,
##    with discovery direction (inner loop, per scenario's med_cl)

pairs_tbl <- read.csv(here("processed-data", "22_Mediation", "03_Mediation_Xenium", "Oligo3_Astro_mediator_outcome_pairs.csv"))
mediation_eval <- read.csv(here("processed-data", "21_Xenium", "06_xenium_ALL_probe_eval", "ECR_mediator_outcome_xenium_eval.csv"))

message(Sys.time(), sprintf(" - %d scenarios to run from pairs table", nrow(pairs_tbl)))


## Outer loop: one scenario (row of pairs_tbl) at a time

all_scenario_results <- pmap(pairs_tbl, function(med_cl, mediator_datatype, outcome_datatype, outcome_cl, med_cl_test) {

    message(Sys.time(), sprintf(" - SCENARIO: med_cl=%s mediator_datatype=%s outcome_datatype=%s outcome_cl=%s med_cl_test=%s",
                                 med_cl, mediator_datatype, outcome_datatype, outcome_cl, med_cl_test))

    sce_pb          <- load_pb_cached(get_pb_fn(outcome_datatype))
    sce_mediator_pb <- load_pb_cached(get_pb_fn(mediator_datatype))

    baseline_tt_out   <- get_baseline_tt_cached(outcome_datatype, outcome_cl)
    dge_base_filtered <- get_filtered_dge_cached(get_pb_fn(outcome_datatype), outcome_cl, baseline_formula_str)

    if (!setequal(rownames(dge_base_filtered), baseline_tt_out$gene_name)) {
        warning(sprintf(
            "Scenario med_cl=%s outcome_cl=%s: re-derived filtered gene set doesn't match saved baseline (%d vs %d genes) -- skipping.",
            med_cl, outcome_cl, nrow(dge_base_filtered), nrow(baseline_tt_out)))
        return(NULL)
    }

    ## candidate mediator genes for THIS scenario's med_cl, restricted to
    ## genes present on this scenario's mediator pb panel
    med_degs <- mediation_eval |>
        filter(med_cl == !!med_cl) |>
        distinct(mediator) |>
        filter(mediator %in% rownames(sce_mediator_pb))

    if (nrow(med_degs) < 1) {
        message(sprintf("  - no candidate mediators for med_cl=%s on this panel -- skipping scenario", med_cl))
        return(NULL)
    }
    message(sprintf("  - %d candidate mediator genes for med_cl=%s", nrow(med_degs), med_cl))

    ## ------------------------------------------------------------------
    ## Inner loop: one candidate mediator gene at a time, within this
    ## scenario. Note med_cl_test (not med_cl) is what's used to pull
    ## expression -- for APOE-stratified scenarios these differ
    ## (e.g. med_cl="Astro.1" but med_cl_test="Astro.1_APOE_high").
    ## ------------------------------------------------------------------
    scenario_results <- map(med_degs$mediator, possibly(function(gene_name) {

        med_vec_full <- get_mediator_vector(sce_mediator_pb, med_cl_test, gene_name)

        out_brnum <- dge_base_filtered$samples$BrNum
        common_donors <- intersect(out_brnum, names(med_vec_full))
        if (length(common_donors) < 10) {
            message(sprintf("  - %s|%s: only %d donors overlap -- skipping", med_cl_test, gene_name, length(common_donors)))
            return(NULL)
        }

        dge_sub <- dge_base_filtered[, out_brnum %in% common_donors]
        dge_sub <- dge_sub[, match(common_donors, out_brnum)]
        med_vec <- med_vec_full[dge_sub$samples$BrNum]

        fit <- fit_mediation_de(dge_sub, baseline_formula_str, batch, med_vec)

        mediation_validation_details(mediation_eval, med_cl, gene_name, baseline_tt_out, fit$tt, fit$ttM, Pvalthr = Pvalthr) |>
            mutate(med_cl_test = med_cl_test, mediator_datatype = mediator_datatype,
                   outcome_datatype = outcome_datatype, outcome_cl = outcome_cl,
                   n_donors = length(common_donors), .before = 1)
    }, otherwise = NULL))

    bind_rows(scenario_results)
})

mediation_summary <- bind_rows(all_scenario_results)

write.csv(mediation_summary,
          file = here(out_dir, sprintf("mediation_summary-all_scenarios_Pval%.2f.csv", Pvalthr)),
          row.names = FALSE)
saveRDS(mediation_summary, file = here(out_dir, sprintf("mediation_summary-all_scenarios_Pval%.2f.rds", Pvalthr)))

message(Sys.time(), sprintf(" - Mediation validation complete: %d scenario x pair rows, %d with mediated=TRUE",
                             nrow(mediation_summary), sum(mediation_summary$mediated, na.rm = TRUE)))

# slurmjobs::job_single('03_Mediation_Xenium_pairs', create_shell = TRUE, memory = '15G', command = "Rscript 03_Mediation_Xenium_pairs.R")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
