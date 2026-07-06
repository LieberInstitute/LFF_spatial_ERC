## Louise Huuki-Myers June 2026 (base voomLmFit script)
## Mediation step added to validate the snRNA-seq mediation screen
## (01_mediation_screening.R) in Xenium spatial pseudobulk data.
##
## Updated to match 06_Clusterwise_voomLmFit_Xenium.R, which now:
##   - takes --datatype on the CLI, selecting between 4 pseudobulk objects
##   - critically, "Xenium_Oligo.3_Astro" / "Xenium_Oligo.3_Astro_SpX" are
##     pseudobulk objects ALREADY RESTRICTED to just Oligo.3 (outcome) and
##     Astro (mediator) clusters -- this is the natural input for mediation,
##     since the outcome cluster and every candidate mediator cluster live
##     in the same sce_pb object with no extra loading needed.
##
## This script does NOT refit the clusterwise baseline DE loop -- it
## assumes 06_Clusterwise_voomLmFit_Xenium.R has already been run
## (unmodified) for the chosen --datatype, and reads its saved per-cluster
## topTable output for the outcome cluster directly from DE_data_dir. It then
## re-derives just the filtered DGEList for that one cluster and fits the
## mediation model (Steps 3a/3b) against a table of candidate mediators.

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

## CLI: mirrors 06_Clusterwise_voomLmFit_Xenium.R's --datatype option so
## this script can be pointed at whichever pseudobulk object you ran that
## script against.

mopt <- matrix(
    c("datatype",   "d", "1", "character", "Data type (must match 06_Clusterwise_voomLmFit_Xenium.R run)",
      "outcome",     "o", "1", "character", "Outcome cluster (default: Oligo.3)",
      "outcome",     "o", "1", "character", "Outcome cluster (default: Oligo.3)",
      "mediatorfdr", "m", "1", "double",    "FDR threshold for calling a gene mediated (default: 0.05)"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(mopt)

## example
## opt$datatype = "cell_type_anno"

if (is.null(opt$datatype)) opt$datatype <- "Xenium_cell_type_anno"
if (is.null(opt$outcome_cluster))  opt$outcome_cluster  <- "Oligo.3"
if (is.null(opt$mediatorfdr) || is.na(opt$mediatorfdr)) opt$mediatorfdr <- 0.05

out_clus <- opt$outcome_cluster
Pvalthr <- opt$mediatorfdr


## pb_fn selection -- kept identical to 06_Clusterwise_voomLmFit_Xenium.R
## so the two scripts never drift apart on which file a given --datatype
## points to.

if (opt$datatype == "Xenium_cell_type_anno") {
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno.RDS")
} else if (opt$datatype == "Xenium_cell_type_anno_SpX") {
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-cell_type_anno_SpX.RDS")
} else if (opt$datatype == "Xenium_Oligo.3_Astro") {
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-Oligo.3_Astro.RDS")
} else if (opt$datatype == "Xenium_Oligo.3_Astro_SpX") {
    pb_fn <- here("processed-data", "21_Xenium", "19_xenium_pseudobulk_DE_prep", "spe_xenium_pseudo_DGE-Oligo.3_Astro_SpX.RDS")
} else {
    stop("non-valid datatype")
}

batch <- "chip"
message(Sys.time(), sprintf(" - Datatype = %s, loading '%s'", opt$datatype, pb_fn))

#### Set up dirs (matches 06_Clusterwise_voomLmFit_Xenium.R's DE_data_dir) ####
DE_data_dir <- here("processed-data", "12_voomLmFit", "06_Clusterwise_voomLmFit_Xenium", sprintf("vlmf_%s", opt$datatype))

## TODO fix output and location of this file
# mediation_dir <- here("processed-data", "22_Mediation", sprintf("vlmf_%s_mediation", opt$datatype))
# if (!dir.exists(mediation_dir)) dir.create(mediation_dir, recursive = TRUE)

# mediation_dir <- here("processed-data", "22_Mediation", sprintf("vlmf_%s_mediation", opt$datatype))
# if (!dir.exists(mediation_dir)) dir.create(mediation_dir, recursive = TRUE)

#### Load the pseudobulk data ####
sce_pb <- readRDS(pb_fn)

dim(sce_pb)
table(sce_pb$registration_variable)

sce_pb$APOE_carrier_syn <- gsub("\\+", "", sce_pb$APOE_carrier)
table(sce_pb$APOE_carrier_syn)

#### PART 1: load the already-computed baseline (Step 1, X->Y) output ####
## from the unmodified 06_Clusterwise_voomLmFit_Xenium.R run for this datatype,
## instead of refitting it here.

baseline_formula_str <- "APOE_carrier_syn + Age + Anc_Afr" ## must match 06_Clusterwise_voomLmFit_Xenium.R exactly

message(Sys.time(), sprintf(" - Mediation validation: datatype = %s, outcome cluster = %s", opt$datatype, out_clus))

baseline_tt_fn <- here(DE_data_dir, sprintf("voomLmFit_%s_%s.rds", opt$datatype, out_clus))

if (!file.exists(baseline_tt_fn)) {
    stop(sprintf("Baseline output not found: %s -- run 06_Clusterwise_voomLmFit_Xenium.R --datatype %s first.",
                 baseline_tt_fn, opt$datatype))
}
baseline_tt_out <- readRDS(baseline_tt_fn)

## Check for a usable gene-id column 
gene_id_col <- "gene_id"
gene_id_col %in% colnames(baseline_tt_out)

## Re-derive the filtered DGEList for the outcome cluster only (the base
## script doesn't save this, only the topTable). One cluster's worth of
## calcNormFactors/filterByExpr -- cheap, and deterministic given the same
## pb_fn + formula, so it reproduces exactly the gene/donor set behind
## baseline_tt_out.

get_filtered_dge <- function(sce_pb, clus, formula_str) {
    dge <- sce_pb[, sce_pb$registration_variable == clus]
    des <- model.matrix(as.formula(paste("~", formula_str)), data = colData(dge))
    dge <- edgeR::calcNormFactors(dge)
    keep <- edgeR::filterByExpr.DGEList(dge, design = as.data.frame(des))
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    edgeR::calcNormFactors(dge)
}
dge_base_filtered <- get_filtered_dge(sce_pb, out_clus, baseline_formula_str)

if (!setequal(rownames(dge_base_filtered), baseline_tt_out$gene_name)) {
    stop(sprintf(
        "Re-derived filtered gene set for %s does not match the saved baseline output (%d vs %d genes). ",
        out_clus, nrow(dge_base_filtered), nrow(baseline_tt_out)),
        "Check that pb_fn, batch, and baseline_formula_str exactly match what produced the saved RDS.")
}

#### PART 2: mediation step (Steps 3a/3b) ####

## Candidate mediators: load from the snRNA-seq mediation screen output
## Xenium is a targeted panel -- restrict candidate mediators to genes actually present on this pseudobulk object's panel.

mediation_eval<- read.csv(here("processed-data", "21_Xenium", "06_xenium_ALL_probe_eval", "ECR_mediator_outcome_xenium_eval.csv"))
    
med_degs <- mediation_eval |>
    select(med_cl, mediator) |>
    unique() |>
    filter(mediator %in% rownames(sce_pb))

message(Sys.time(), sprintf(" - %d candidate mediators present on Xenium panel", nrow(med_degs)))
if (nrow(med_degs) < 1) stop("No candidate mediator genes from the snRNA-seq screen are present on the Xenium panel.")

## TODO remove ?
# ## If using an _SpX datatype, med_degs$cluster (plain names like "Astro.2")
# ## may not match sce_pb$registration_variable (possibly domain-crossed
# ## labels) -- see the SpX NOTE at the top of this file before proceeding.
# if (grepl("_SpX$", opt$datatype)) {
#     unmatched_clusters <- setdiff(unique(med_degs$cluster), unique(as.character(sce_pb$registration_variable)))
#     if (length(unmatched_clusters) > 0) {
#         warning(sprintf(
#             "datatype '%s' is a SpX (domain-crossed) object, and %d mediator cluster label(s) don't exactly match sce_pb$registration_variable: %s. ",
#             opt$datatype, length(unmatched_clusters), paste(head(unmatched_clusters, 5), collapse = ", ")),
#             "See the SpX NOTE at the top of this script -- get_mediator_vector() will silently return NULL for these and they'll be skipped.")
#     }
# }


## Build a BrNum-indexed mediator expression vector from the pseudobulk
## object itself, for one (cluster, gene_id) pair at a time. Since the
## Xenium_Oligo.3_Astro(_SpX) objects already contain both the outcome
## and mediator clusters, this is a straight subset of sce_pb -- no
## separate/larger object needs to be loaded.

get_mediator_vector <- function(sce_pb, med_cluster, gene_id) {
    
    stopifnot(gene_id %in% rownames(sce_pb))
    clus_sce <- sce_pb[gene_id, sce_pb$registration_variable == med_cluster]
    brnum <- as.character(colData(clus_sce)$BrNum)
    
    stopifnot(!anyDuplicated(brnum))

    mediator_vector <- as.numeric(assay(clus_sce, "logcounts"))
    names(mediator_vector) <- brnum
    return(mediator_vector)
}

# med_vec_test <- get_mediator_vector(sce_pb, med_cluster = "Astro.2", gene_id = "FZD8")

mediation_eval |> filter(mediator == "FZD8")

#### Core mediation fit ####
## same voomLmFit -> eBayes -> topTable pipeline as
## the base script, with the mediator added to the design. The design
## keeps its intercept, so the carrier coefficient keeps its usual name
## ("APOE_carrier_E4") -- no makeContrasts needed.

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

fit_mediation_de_test <- fit_mediation_de(dge_base_filtered, baseline_formula_str, batch_col = batch,
                                          med_vec = get_mediator_vector(sce_pb, med_cluster = "Astro.1", gene_id = "PLPPR4"))

                                                                                            ## gene_nameP.Value   t
# tt_baseline |> filter(gene_name == "TESPA1") |> select(gene_name, P.Value, t)               ## TESPA1     0.0618  1.99  # validates in baseline
# fit_mediation_de_test$tt |> filter(gene_name == "TESPA1") |> select(gene_name, P.Value, t)  ## TESPA1     0.0961  1.77  # controlling for "FZD8" lessens effect but still signif
# fit_mediation_de_test$ttM |> filter(gene_name == "TESPA1") |> select(gene_name, P.Value, t) ## TESPA1      0.918  0.11  # not associated with "FZD8" expression

## Step 3b classification: attenuated-and-M-associated = "mediated"

classify_mediation <- function(tt_baseline, tt_mediation, ttM_mediation, Pvalthr = 0.10) {
    base_sig    <- tt_baseline  |> filter(P.Value < Pvalthr) |> pull(gene_name)
    carrier_sig <- tt_mediation |> filter(P.Value < Pvalthr) |> pull(gene_name)
    med_sig     <- ttM_mediation |> filter(P.Value < Pvalthr) |> pull(gene_name)

    lost     <- setdiff(base_sig, carrier_sig)
    gained   <- setdiff(carrier_sig, base_sig)
    mediated <- intersect(lost, med_sig)

    tibble(
        baseline_n = length(base_sig), 
        carrier_n = length(carrier_sig), 
        med_n = length(med_sig),
        lost_n = length(lost), 
        gained_n = length(gained), 
        mediated_n = length(mediated),
        pass_both = mediated_n > 0,
        mediated_genes = list(mediated)
    )
}

test_classify_mediation <- classify_mediation(tt_baseline = baseline_tt_out,
                                               tt_mediation = fit_mediation_de_test$tt,
                                               ttM_mediation = fit_mediation_de_test$ttM)


mediation_validation_details <- function(med_cluster, 
                                         gene_name, 
                                         tt_baseline, 
                                         tt_mediation, 
                                         ttM_mediation,
                                         Pvalthr = 0.10){
    
    mediation_eval |>
        filter(med_cl == med_cluster, mediator == gene_name) |>
        select(med_cl:t_med)|>
        rename_with(~paste0(.x, "_SN"), fdr:t_med) |>
        left_join(tt_baseline |> select(outcome = gene_name, P.Value_base = P.Value, t_base = t), by = join_by(outcome)) |>
        left_join(tt_mediation |> select(outcome = gene_name, P.Value_carrier = P.Value, t_carrier = t), by = join_by(outcome)) |>
        left_join(ttM_mediation |> select(outcome = gene_name, P.Value_med = P.Value, t_med = t), by = join_by(outcome)) |>
        mutate(base_sig = P.Value_base < Pvalthr,
               base_dir = sign(t_SN) == sign(t_base),
               base_valid = base_sig & base_dir,
               carrier_sig = P.Value_carrier < Pvalthr,
               carrier_dir = sign(t_carrier) == sign(t_base),
               med_sig = P.Value_med < Pvalthr,
               mediated = !carrier_sig & med_sig)
    
}

baseline_tt_out |> filter(gene_name %in% 
                              (mediation_eval |> filter(med_cl == "Astro.1", mediator == "PLPPR4") |> pull(outcome))
                          ) 

## Loop over candidate mediator genes, fit mediation model, classify

message(Sys.time(), sprintf(" - Running mediation model for %d candidate mediators (Pvalthr=%.4f)", nrow(med_degs), Pvalthr))

mediation_results <- pmap(
    list(med_degs$med_cl, med_degs$mediator),
    possibly(function(med_cluster, gene_name) {

        med_vec_full <- get_mediator_vector(sce_pb, med_cluster, gene_name)
        if (is.null(med_vec_full)) return(NULL)

        out_brnum <- dge_base_filtered$samples$BrNum
        common_donors <- intersect(out_brnum, names(med_vec_full))
        if (length(common_donors) < 10) return(NULL)

        dge_sub <- dge_base_filtered[, out_brnum %in% common_donors]
        dge_sub <- dge_sub[, match(common_donors, out_brnum)]
        med_vec <- med_vec_full[out_brnum]

        fit <- fit_mediation_de(dge_sub, baseline_formula_str, batch, med_vec)

        run_label <- sprintf("%s|%s", med_cluster, gene_name)
        # saveRDS(fit, file = here(mediation_dir, sprintf("mediation_%s_%s_%s.rds",
        #                                                  gsub("\\.", "", med_cluster), gene_id, gene_name)))

        summary_row <- mediation_validation_details(med_cluster, gene_name, baseline_tt_out, fit$tt, fit$ttM, Pvalthr = Pvalthr)
        # summary_row |> mutate(run = run_label, med_cluster = med_cluster,
        #                        gene_name = gene_name, n_donors = length(common_donors), .before = 1)
        
        return(summary_row)
    }, otherwise = NULL)
)

mediation_summary <- bind_rows(mediation_results)

write.csv(mediation_summary |> select(-mediated_genes),
          file = here(mediation_dir, sprintf("mediation_summary-%s_%s.csv", opt$datatype, out_clus)),
          row.names = FALSE)

# saveRDS(mediation_summary, file = here(mediation_dir, sprintf("mediation_summary-%s_%s.rds", opt$datatype, out_clus)))

message(Sys.time(), sprintf(" - Mediation validation complete: %d/%d mediators tested, %d with pass_both=TRUE",
                             nrow(mediation_summary), nrow(med_degs), sum(mediation_summary$pass_both, na.rm = TRUE)))

# slurmjobs::job_single('07_Mediation_Xenium', create_shell = TRUE, memory = '25G',
#                       command = "Rscript 07_Mediation_Xenium.R --datatype Xenium_Oligo.3_Astro")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
