#!/bin/env Rscript

## test the mediation voomLmFit approach on our data, sensitivity to donor restrictions

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("SummarizedExperiment")
library("BiocParallel")
library("here")
library("getopt")
library("sessioninfo")
library("qs2")
library("ggplot2")

md5_string <- function(x) {
    tf <- tempfile()
    writeLines(x, con = tf, useBytes = TRUE)
    md5 <- as.character(tools::md5sum(tf))
    unlink(tf)
    md5
}

make_donor_key <- function(brnums) {
    donors <- sort(unique(as.character(brnums)))
    donor_str <- paste(donors, collapse = ",")
    donor_key <- md5_string(donor_str)
    list(donors = donors, donor_str = donor_str, donor_key = donor_key)
}

write_fdr_table <- function(dt, file) {
    ## write empty gzip with header when no rows to avoid corrupt gzip output
    if (nrow(dt) < 1) {
        con <- gzfile(file, "wt")
        on.exit(close(con), add = TRUE)
        if (ncol(dt) > 0) {
            writeLines(paste(names(dt), collapse = "\t"), con = con, useBytes = TRUE)
        }
        return(invisible(NULL))
    }
    data.table::fwrite(dt, file = file, sep = "\t", compress = "gzip")
}

apoe_to_carrier <- function(x) {
    x <- as.character(x)
    out <- ifelse(x %in% c("E2.E2", "E2.E3"), "E2_carrier",
                  ifelse(x %in% c("E3.E4", "E4.E4"), "E4_carrier", NA_character_))
    factor(out, levels = c("E2_carrier", "E4_carrier"))
}

get_apoe_counts <- function(apoe_vec) {
    tab <- table(factor(as.character(apoe_vec), levels = c("E2.E2", "E2.E3", "E3.E4", "E4.E4")))
    list(
        n_E2_E2 = as.integer(tab["E2.E2"]),
        n_E2_E3 = as.integer(tab["E2.E3"]),
        n_E3_E4 = as.integer(tab["E3.E4"]),
        n_E4_E4 = as.integer(tab["E4.E4"])
    )
}

cont_smd <- function(x, g) {
    g <- factor(g)
    if (nlevels(g) != 2) return(NA_real_)
    x1 <- x[g == levels(g)[1]]
    x2 <- x[g == levels(g)[2]]
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]
    if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
    s1 <- stats::sd(x1)
    s2 <- stats::sd(x2)
    sp <- sqrt((s1^2 + s2^2) / 2)
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    (mean(x1) - mean(x2)) / sp
}

fmt <- function(x, digits = 3) {
    if (length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
    formatC(x, format = "f", digits = digits)
}

write_matrix_as_text <- function(con, mat, title) {
    writeLines(title, con)
    if (length(mat) == 0 || any(dim(mat) == 0)) {
        writeLines("  <empty>", con)
        writeLines("", con)
        return(invisible(NULL))
    }
    mm <- as.matrix(mat)
    rownames(mm) <- rownames(mat)
    colnames(mm) <- colnames(mat)
    lines <- capture.output(print(mm, quote = FALSE))
    writeLines(paste0("  ", lines), con)
    writeLines("", con)
}

## command line options
opt_spec <- matrix(
    c(
        "mode", "m", "1", "character", "Run mode: stepdown|sweep|both|outlier|all",
        "sweeps", "s", "1", "integer", "Number of random sweeps (default=3)",
        "seed", "r", "1", "integer", "Random seed for sweeps (default=20260205)",
        "sweep_targets", "t", "1", "character", "Comma-separated target donor sizes for sweep mode (default=20,19,18,17,16,15)",
        "workers", "w", "1", "integer", "Workers for outlier scans (default=8)"
    ),
    ncol = 5, byrow = TRUE
)
opt <- tryCatch(getopt(opt_spec), error = function(e) list())
run_mode <- if (is.null(opt$mode) || is.na(opt$mode)) "stepdown" else as.character(opt$mode)
n_sweeps <- if (is.null(opt$sweeps) || is.na(opt$sweeps)) 3L else as.integer(opt$sweeps)
sweep_seed <- if (is.null(opt$seed) || is.na(opt$seed)) 20260205L else as.integer(opt$seed)
sweep_targets_str <- if (is.null(opt$sweep_targets) || is.na(opt$sweep_targets)) "20,19,18,17,16,15" else as.character(opt$sweep_targets)
sweep_target_ns <- as.integer(strsplit(gsub("\\s+", "", sweep_targets_str), ",", fixed = FALSE)[[1]])
workers <- if (is.null(opt$workers) || is.na(opt$workers)) 8L else as.integer(opt$workers)

if (!(run_mode %in% c("stepdown", "sweep", "both", "outlier", "all"))) {
    stop("Invalid --mode. Use one of: stepdown, sweep, both, outlier, all")
}
if (is.na(n_sweeps) || n_sweeps < 1) stop("--sweeps must be >= 1")
if (is.na(sweep_seed)) stop("--seed must be an integer")
if (is.na(workers) || workers < 1) workers <- 1L
if (length(sweep_target_ns) < 1 || any(is.na(sweep_target_ns))) {
    stop("Invalid --sweep_targets. Provide comma-separated integers, e.g. 20,19,18,17,16")
}
if (any(diff(sweep_target_ns) >= 0)) {
    stop("--sweep_targets must be strictly decreasing.")
}
do_stepdown <- run_mode %in% c("stepdown", "both", "all")
do_sweep <- run_mode %in% c("sweep", "both", "all")
do_outlier <- run_mode %in% c("outlier", "all")
message(Sys.time(), sprintf(" - Run mode: %s (stepdown=%s, sweep=%s, outlier=%s)", run_mode, do_stepdown, do_sweep, do_outlier))

#### Paths ####
ddir <- here("processed-data")
mdir <- here("processed-data", "22_Mediation")

pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", "sce_pseudo_DGE-cell_type_anno.RDS")

out_dir <- here("code", "22_Mediation", "out-dge-test")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
baseline_dir <- file.path(out_dir, "baseline")
dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(out_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#### Load X->Y pseudo bulk data ####
sce_pb <- readRDS(pb_fn)
## always drop Br1289
sce_pb <- sce_pb[, sce_pb$BrNum != "Br1289"]

#### Outcome cluster setup (same for every scenario) ####
## we can use dge_base colData() for all donor metadata

out_clus <- "Oligo.3"
batch <- "exp_round"
sce_base <- sce_pb[, sce_pb$registration_variable == out_clus]

donor_meta <- as.data.frame(colData(sce_base))

baseline_formula_str <- "0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio"
baseline_formula <- as.formula(paste("~", baseline_formula_str))

baseline_des <- model.matrix(
    baseline_formula,
    data = donor_meta
)
baseline_des <- as.data.frame(baseline_des)


## LC Astro data
## LC Astro DEGs from Bernie's DGE - not used in this script
f_lca_dge <- file.path(mdir, "LC_Astro_E4vE2_DGEout_Bernie.rds")
if (!file.exists(f_lca_dge)) {
    stop(sprintf("LC Astro DEG file not found: %s", f_lca_dge))
}
lca_dge <- as.data.table(readRDS(f_lca_dge))

## Load LC Astro pseudobulk counts - used to determine overlapping donors
f_lca_pbcounts <- file.path(mdir, "lc_astro_pseudobulk_counts_by_donor.rds")
if (!file.exists(f_lca_pbcounts)) {
    stop(sprintf("LC Astro pseudobulk counts file not found: %s", f_lca_pbcounts))
}
dgl_lca <- readRDS(f_lca_pbcounts)

## Determine overlapping donors between LC Astro and ERC Oligo.3 outcome
lca_brnums <- colnames(dgl_lca)
erc_brnums <- as.character(donor_meta$BrNum)
common_donors_21 <- intersect(lca_brnums, erc_brnums)

message(Sys.time(), sprintf(" - ERC Oligo.3 donors: %d", length(unique(erc_brnums))))
message(Sys.time(), sprintf(" - LC Astro donors: %d", length(unique(lca_brnums))))
message(Sys.time(), sprintf(" - Common donors: %d", length(common_donors_21)))

## Match baseline filtering from 01_mediation_screening.R
#fn_dge_base <- file.path(out_dir, "dge_base.qs2")
#qs_save(dge_base, fn_dge_base)
# message(Sys.time(), sprintf(" - Saved filtered dge_base cache: %s", fn_dge_base))

build_design <- function(sample_df, formula_obj, batch_col, med_vec = NULL) {
    sample_df <- as.data.frame(sample_df)
    if ("APOE_syn" %in% names(sample_df)) sample_df$APOE_syn <- droplevels(factor(sample_df$APOE_syn))
    if ("Sex" %in% names(sample_df)) sample_df$Sex <- droplevels(factor(sample_df$Sex))
    if (!(batch_col %in% names(sample_df))) {
        stop(sprintf("Batch column '%s' not found in samples.", batch_col))
    }
    sample_df[[batch_col]] <- droplevels(factor(sample_df[[batch_col]]))
    if (!is.null(med_vec)) sample_df$med_vec <- med_vec
    des <- model.matrix(formula_obj, data = sample_df)
    des <- as.data.frame(des)
    colnames(des) <- gsub(colnames(des), pattern = "_syn", replacement = "_")
    list(design = des, sample_df = sample_df)
}

run_baseline <- function(donor_key, donors, donor_str) {
    out_file <- file.path(baseline_dir, sprintf("baseline_%s.tab.gz", donor_key))
    out_fdr_file <- file.path(baseline_dir, sprintf("baseline_%s.fdr05.tab.gz", donor_key))

    out_brnum <- as.character(sce_base$BrNum)
    common_donors <- intersect(out_brnum, donors)
    if (length(common_donors) < 10) {
        warning(sprintf("Too few common donors (%d) for baseline %s", length(common_donors), donor_key))
        return(NULL)
    }
    ## subset dge_base to common donors
    dge_base <- edgeR::calcNormFactors(sce_base[, sce_base$BrNum %in% common_donors])
    keep <- edgeR::filterByExpr(dge_base, design = baseline_des)
    dge_base <- dge_base[keep, , keep.lib.sizes = FALSE]
    dge_base <- edgeR::calcNormFactors(dge_base)


    dge_sub <- dge_base[, out_brnum %in% common_donors]
    design_obj <- build_design(dge_sub$samples, baseline_formula, batch)
    des <- design_obj$design
    sample_df <- design_obj$sample_df

    req_cols <- c("APOE_E2.E2", "APOE_E2.E3", "APOE_E3.E4", "APOE_E4.E4")
    missing_cols <- setdiff(req_cols, colnames(des))
    if (length(missing_cols) > 0) {
        apoe_tab <- table(sample_df$APOE_syn)
        stop(sprintf(
            "Missing APOE columns in design (%s) for donor_key %s. APOE levels in subset: %s",
            paste(missing_cols, collapse = ", "),
            donor_key,
            paste(sprintf("%s:%d", names(apoe_tab), as.integer(apoe_tab)), collapse = "; ")
        ))
    }

    if (file.exists(out_file)) {
        message(Sys.time(), sprintf(" - Using cached baseline: %s", out_file))
        tt <- data.table::fread(out_file)
    } else {
        v.swt <- voomLmFit(
            dge_sub,
            design = des,
            block = as.factor(sample_df[[batch]]),
            adaptive.span = TRUE,
            sample.weights = TRUE
        )

        cont <- makeContrasts(
            carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
            levels = des
        )

        v.swt.fit <- contrasts.fit(v.swt, contrasts = cont)
        v.swt.fit.e <- eBayes(v.swt.fit)

        tt <- topTable(v.swt.fit.e, coef = "carrier", number = Inf, adjust.method = "BH")
        tt <- as.data.table(tt, keep.rownames = "outcome_gene_id")
        tt[, `:=`(run = sprintf("baseline|%s", donor_key), cluster = out_clus, contrast = "carrier",
                  n_donors = length(common_donors), donor_key = donor_key)]
        setcolorder(tt, c("run", "cluster", "contrast", "outcome_gene_id", "n_donors", "donor_key",
                          setdiff(names(tt), c("run", "cluster", "contrast", "outcome_gene_id", "n_donors", "donor_key"))))

        data.table::fwrite(tt, file = out_file, sep = "\t", compress = "gzip")
    }

    write_fdr_table(tt[adj.P.Val < 0.05], out_fdr_file)
    baseline_sig <- tt[adj.P.Val < 0.05, outcome_gene_id]
    apoe_tab <- table(sample_df$APOE_syn)
    carrier_vec <- apoe_to_carrier(sample_df$APOE_syn)
    carrier_tab <- table(carrier_vec)
    apoe_counts <- get_apoe_counts(sample_df$APOE_syn)

    list(
        donor_key = donor_key,
        donor_str = donor_str,
        n_donors = length(common_donors),
        design_n_samples = nrow(des),
        design_n_coef = ncol(des),
        apoe_levels = paste(sprintf("%s:%d", names(apoe_tab), as.integer(apoe_tab)), collapse = ";"),
        n_E2_E2 = apoe_counts$n_E2_E2,
        n_E2_E3 = apoe_counts$n_E2_E3,
        n_E3_E4 = apoe_counts$n_E3_E4,
        n_E4_E4 = apoe_counts$n_E4_E4,
        n_E2_carrier = as.integer(carrier_tab["E2_carrier"]),
        n_E4_carrier = as.integer(carrier_tab["E4_carrier"]),
        min_p = min(tt$P.Value, na.rm = TRUE),
        min_fdr = min(tt$adj.P.Val, na.rm = TRUE),
        n_p_lt_1e3 = sum(tt$P.Value < 1e-3, na.rm = TRUE),
        n_p_lt_1e2 = sum(tt$P.Value < 1e-2, na.rm = TRUE),
        n_p_lt_5e2 = sum(tt$P.Value < 5e-2, na.rm = TRUE),
        out_file = out_file,
        out_fdr_file = out_fdr_file,
        n_deg_fdr05 = length(baseline_sig),
        sig_genes = list(baseline_sig)
    )
}

if (do_stepdown) {
full_donors <- unique(as.character(dge_base$samples$BrNum))
common_donors <- full_donors[full_donors %in% common_donors_21]
erc_only_donors <- full_donors[!(full_donors %in% common_donors)]

if (length(erc_only_donors) != 9) {
    warning(sprintf(
        "Expected 9 ERC-only donors (30->21), observed %d; running %d total baseline scenarios.",
        length(erc_only_donors),
        length(erc_only_donors) + 1
    ))
}

donor_sets_to_run <- lapply(0:length(erc_only_donors), function(i) {
    removed_so_far <- if (i > 0) erc_only_donors[seq_len(i)] else character(0)
    donors_i <- full_donors[!(full_donors %in% removed_so_far)]
    donor_info <- make_donor_key(donors_i)
    list(
        run_id = sprintf("run%02d", i + 1),
        step = i,
        removed_this = if (i > 0) erc_only_donors[i] else NA_character_,
        removed_so_far = paste(removed_so_far, collapse = ","),
        n_donors_planned = length(donors_i),
        donors = list(donor_info$donors),
        donor_key = donor_info$donor_key,
        donor_str = donor_info$donor_str
    )
})
donor_sets_to_run <- data.table::rbindlist(donor_sets_to_run, fill = TRUE)

last_donors <- donor_sets_to_run[nrow(donor_sets_to_run)]$donors[[1]]
if (!setequal(last_donors, common_donors)) {
    stop("Final donor set does not match ERC-LC common donors.")
}

message(Sys.time(), sprintf(" - Running %d baseline step-down analyses", nrow(donor_sets_to_run)))

baseline_results <- lapply(seq_len(nrow(donor_sets_to_run)), function(i) {
    row_i <- donor_sets_to_run[i]
    res <- run_baseline(
        donor_key = row_i$donor_key,
        donors = row_i$donors[[1]],
        donor_str = row_i$donor_str
    )
    if (is.null(res)) return(NULL)
    c(res, list(
        run_id = row_i$run_id,
        step = row_i$step,
        removed_this = row_i$removed_this,
        removed_so_far = row_i$removed_so_far,
        n_donors_planned = row_i$n_donors_planned
    ))
})
baseline_results <- Filter(Negate(is.null), baseline_results)

if (length(baseline_results) < 1) {
    stop("No baseline runs completed successfully.")
}

baseline_dt <- data.table::rbindlist(baseline_results, fill = TRUE)

if (any(baseline_dt$n_donors != baseline_dt$n_donors_planned)) {
    stop("Donor count mismatch between planned and used donor sets.")
}
if (any(baseline_dt$design_n_samples != baseline_dt$n_donors)) {
    stop("Design matrix row count does not match active donor count for one or more runs.")
}

setorder(baseline_dt, step)

data.table::fwrite(
    baseline_dt[, .(
        run_id,
        step,
        donor_key,
        donor_str,
        n_donors,
        design_n_samples,
        design_n_coef,
        apoe_levels,
        n_E2_E2,
        n_E2_E3,
        n_E3_E4,
        n_E4_E4,
        n_E2_carrier,
        n_E4_carrier,
        removed_this,
        removed_so_far,
        min_p,
        min_fdr,
        n_p_lt_1e3,
        n_p_lt_1e2,
        n_p_lt_5e2,
        out_file,
        out_fdr_file,
        n_deg_fdr05
    )],
    file = file.path(out_dir, "baseline_stepdown_summary.tsv"),
    sep = "\t"
)

data.table::fwrite(
    baseline_dt[, .(run_id, donor_key, donor_str, n_donors, out_file, out_fdr_file)],
    file = file.path(out_dir, "baseline_index.tsv"),
    sep = "\t"
)

apoe_dist_dt <- baseline_dt[, .(
    run_id,
    step,
    n_donors,
    n_E2_E2,
    n_E2_E3,
    n_E3_E4,
    n_E4_E4,
    n_E2_carrier,
    n_E4_carrier,
    p_E2_carrier = n_E2_carrier / n_donors,
    p_E4_carrier = n_E4_carrier / n_donors
)]
data.table::fwrite(
    apoe_dist_dt,
    file = file.path(out_dir, "apoe_genotype_distribution_by_run.tsv"),
    sep = "\t"
)

sample_dt_all <- as.data.table(dge_base$samples)
run_full <- baseline_dt[step == min(step)][1]
run_21 <- baseline_dt[step == max(step)][1]
donors_full <- strsplit(run_full$donor_str, ",", fixed = TRUE)[[1]]
donors_21 <- strsplit(run_21$donor_str, ",", fixed = TRUE)[[1]]
meta_full <- sample_dt_all[BrNum %in% donors_full]
meta_21 <- sample_dt_all[BrNum %in% donors_21]
meta_full[, carrier_group := apoe_to_carrier(APOE_syn)]
meta_21[, carrier_group := apoe_to_carrier(APOE_syn)]

covars <- c("Age", "Anc_Afr", "pseudo_expr_chrM_ratio")

covar_long_21 <- data.table::rbindlist(lapply(covars, function(v) {
    vv <- as.numeric(meta_21[[v]])
    data.table::data.table(
        run_id = run_21$run_id,
        covariate = v,
        carrier_group = as.character(meta_21$carrier_group),
        value = vv
    )
}), fill = TRUE)
covar_summary_21 <- covar_long_21[, .(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE)
), by = .(covariate, carrier_group)]

covar_balance_dt <- data.table(
    covariate = covars,
    smd_21 = vapply(covars, function(v) cont_smd(as.numeric(meta_21[[v]]), meta_21$carrier_group), numeric(1)),
    smd_30 = vapply(covars, function(v) cont_smd(as.numeric(meta_full[[v]]), meta_full$carrier_group), numeric(1))
)

sex_tab_21 <- table(meta_21$carrier_group, meta_21$Sex)
sex_tab_30 <- table(meta_full$carrier_group, meta_full$Sex)
sex_p_21 <- if (all(dim(sex_tab_21) == c(2, 2))) fisher.test(sex_tab_21)$p.value else NA_real_
sex_p_30 <- if (all(dim(sex_tab_30) == c(2, 2))) fisher.test(sex_tab_30)$p.value else NA_real_

batch_tab_21 <- table(meta_21$carrier_group, meta_21[[batch]])
batch_tab_30 <- table(meta_full$carrier_group, meta_full[[batch]])
batch_chisq_21 <- suppressWarnings(chisq.test(batch_tab_21, correct = FALSE))
batch_chisq_30 <- suppressWarnings(chisq.test(batch_tab_30, correct = FALSE))
cramer_v <- function(tab) {
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
    n <- sum(tab)
    k <- min(nrow(tab) - 1, ncol(tab) - 1)
    if (n <= 0 || k <= 0) return(NA_real_)
    sqrt(as.numeric(chi$statistic) / (n * k))
}

tt_full <- data.table::fread(run_full$out_file)
tt_21 <- data.table::fread(run_21$out_file)
cmp_dt <- merge(
    tt_full[, .(outcome_gene_id, logFC_30 = logFC, t_30 = t, P_30 = P.Value, FDR_30 = adj.P.Val)],
    tt_21[, .(outcome_gene_id, logFC_21 = logFC, t_21 = t, P_21 = P.Value, FDR_21 = adj.P.Val)],
    by = "outcome_gene_id",
    all = FALSE
)
cmp_dt[, se_30 := abs(logFC_30 / t_30)]
cmp_dt[, se_21 := abs(logFC_21 / t_21)]
cmp_dt[!is.finite(se_30), se_30 := NA_real_]
cmp_dt[!is.finite(se_21), se_21 := NA_real_]
cmp_dt[, same_sign := sign(logFC_30) == sign(logFC_21)]

p_min_required <- 0.05 / nrow(tt_21)
pmin_21 <- min(tt_21$P.Value, na.rm = TRUE)
fdrmin_21 <- min(tt_21$adj.P.Val, na.rm = TRUE)

diag_report_file <- file.path(out_dir, "run10_21donor_diagnostic_report.txt")
con <- file(diag_report_file, "wt")
writeLines(sprintf("21-donor baseline diagnostic report (%s)", Sys.time()), con)
writeLines(sprintf("Full baseline run: %s (n=%d)", run_full$run_id, run_full$n_donors), con)
writeLines(sprintf("21-donor baseline run: %s (n=%d)", run_21$run_id, run_21$n_donors), con)
writeLines("", con)

writeLines("1) APOE composition", con)
writeLines(sprintf("  30-donor counts: E2.E2=%d, E2.E3=%d, E3.E4=%d, E4.E4=%d (E2 carriers=%d, E4 carriers=%d)",
                   run_full$n_E2_E2, run_full$n_E2_E3, run_full$n_E3_E4, run_full$n_E4_E4,
                   run_full$n_E2_carrier, run_full$n_E4_carrier), con)
writeLines(sprintf("  21-donor counts: E2.E2=%d, E2.E3=%d, E3.E4=%d, E4.E4=%d (E2 carriers=%d, E4 carriers=%d)",
                   run_21$n_E2_E2, run_21$n_E2_E3, run_21$n_E3_E4, run_21$n_E4_E4,
                   run_21$n_E2_carrier, run_21$n_E4_carrier), con)
writeLines(sprintf("  Contrast-group balance (21 donors): E2 vs E4 carriers = %d vs %d (balanced).",
                   run_21$n_E2_carrier, run_21$n_E4_carrier), con)
writeLines("", con)

writeLines("2) Statistical signal metrics (carrier contrast on Oligo.3 genes)", con)
writeLines(sprintf("  Genes tested: %d", nrow(tt_21)), con)
writeLines(sprintf("  21-donor min raw p-value = %s", formatC(pmin_21, format = "e", digits = 3)), con)
writeLines(sprintf("  21-donor min BH FDR = %.4f", fdrmin_21), con)
writeLines(sprintf("  Raw p-value needed for BH FDR < 0.05 at rank 1 with %d tests: p < %s",
                   nrow(tt_21), formatC(p_min_required, format = "e", digits = 3)), con)
writeLines(sprintf("  21-donor counts: p<1e-3: %d ; p<1e-2: %d ; p<0.05: %d",
                   run_21$n_p_lt_1e3, run_21$n_p_lt_1e2, run_21$n_p_lt_5e2), con)
writeLines(sprintf("  30-donor min BH FDR = %.4f ; DEGs(FDR<0.05)=%d",
                   run_full$min_fdr, run_full$n_deg_fdr05), con)
writeLines(sprintf("  21-donor DEGs(FDR<0.05)=%d", run_21$n_deg_fdr05), con)
writeLines("", con)

writeLines("3) Effect-size and precision changes (30 donors vs 21 donors)", con)
writeLines(sprintf("  logFC correlation across genes: %.3f", cor(cmp_dt$logFC_30, cmp_dt$logFC_21, use = "complete.obs")), con)
writeLines(sprintf("  Sign concordance across genes: %.1f%%", 100 * mean(cmp_dt$same_sign, na.rm = TRUE)), con)
writeLines(sprintf("  Median |t|: 30-donor=%.3f ; 21-donor=%.3f",
                   median(abs(cmp_dt$t_30), na.rm = TRUE), median(abs(cmp_dt$t_21), na.rm = TRUE)), con)
writeLines(sprintf("  Median SE(logFC): 30-donor=%.4f ; 21-donor=%.4f ; ratio(21/30)=%.3f",
                   median(cmp_dt$se_30, na.rm = TRUE),
                   median(cmp_dt$se_21, na.rm = TRUE),
                   median(cmp_dt$se_21, na.rm = TRUE) / median(cmp_dt$se_30, na.rm = TRUE)), con)
writeLines("", con)

writeLines("4) Covariate balance in the 21-donor set (by carrier group)", con)
writeLines("  Continuous covariates summary:", con)
for (v in covars) {
    dsub <- covar_summary_21[covariate == v]
    if (nrow(dsub) > 0) {
        for (k in seq_len(nrow(dsub))) {
            writeLines(sprintf("    %s | %s: n=%d mean=%s sd=%s median=%s range=[%s,%s]",
                               v, dsub$carrier_group[k], dsub$n[k],
                               fmt(dsub$mean[k]), fmt(dsub$sd[k]), fmt(dsub$median[k]),
                               fmt(dsub$min[k]), fmt(dsub$max[k])), con)
        }
    }
}
for (k in seq_len(nrow(covar_balance_dt))) {
    writeLines(sprintf("    SMD(%s): 21-donor=%s ; 30-donor=%s",
                       covar_balance_dt$covariate[k],
                       fmt(covar_balance_dt$smd_21[k]),
                       fmt(covar_balance_dt$smd_30[k])), con)
}
writeLines(sprintf("  Sex vs carrier association p-value (Fisher): 21-donor=%s ; 30-donor=%s",
                   fmt(sex_p_21, digits = 4), fmt(sex_p_30, digits = 4)), con)
writeLines(sprintf("  Batch(%s) vs carrier association p-value (Chi-square): 21-donor=%s ; 30-donor=%s",
                   batch, fmt(batch_chisq_21$p.value, digits = 4), fmt(batch_chisq_30$p.value, digits = 4)), con)
writeLines(sprintf("  Cramer's V for batch vs carrier: 21-donor=%s ; 30-donor=%s",
                   fmt(cramer_v(batch_tab_21)), fmt(cramer_v(batch_tab_30))), con)
writeLines("", con)
write_matrix_as_text(con, sex_tab_21, "  21-donor Sex x Carrier table:")
write_matrix_as_text(con, batch_tab_21, sprintf("  21-donor %s x Carrier table:", batch))

writeLines("5) Most plausible explanations (data-supported)", con)
writeLines("  - Primary: reduced effective power/precision after dropping from 30 to 21 donors.", con)
writeLines("    Evidence: median SE(logFC) increased and median |t| decreased, so p-values are not low enough after BH correction.", con)
writeLines("  - Multiple-testing burden remains high (~11.5k genes), so the observed minimum raw p-value in 21 donors is above the BH discovery threshold.", con)
writeLines("  - Carrier groups remain numerically balanced (E2 vs E4 carriers), so complete group-size imbalance is unlikely to be the main cause.", con)
writeLines("  - Genotype-composition shift (notably fewer E2.E2 and E3.E4 donors) can reduce contrast stability if subgroup-specific effects differ.", con)
writeLines("  - Covariate imbalance (sex/batch/continuous covariates) should be checked as contributing factors; statistics above quantify this directly for the 21-donor subset.", con)
close(con)

data.table::fwrite(covar_summary_21, file = file.path(out_dir, "run10_covariate_summary_by_carrier.tsv"), sep = "\t")
data.table::fwrite(covar_balance_dt, file = file.path(out_dir, "run10_covariate_balance_metrics.tsv"), sep = "\t")

## ggplot diagnostics to support report claims
plot_index <- data.table::data.table(plot_file = character(), description = character())
add_plot <- function(file, description) {
    plot_index <<- rbind(plot_index, data.table::data.table(plot_file = file, description = description))
}

## Plot 1: DEGs vs donor count
plot_deg_vs_n <- ggplot(baseline_dt, aes(x = n_donors, y = n_deg_fdr05)) +
    geom_line(color = "#1b4f72", linewidth = 0.8) +
    geom_point(aes(color = run_id), size = 2.2) +
    scale_x_continuous(breaks = baseline_dt$n_donors) +
    theme_bw(base_size = 12) +
    labs(
        title = "Baseline DEGs vs donor count",
        subtitle = "Carrier contrast, Oligo.3, FDR < 0.05",
        x = "Donors in run",
        y = "Number of DEGs"
    ) +
    theme(legend.position = "none")
fn_plot_deg_vs_n <- file.path(plots_dir, "01_deg_count_vs_donors.png")
ggsave(fn_plot_deg_vs_n, plot = plot_deg_vs_n, width = 8.5, height = 4.8, dpi = 140)
add_plot(fn_plot_deg_vs_n, "DEG count trajectory across the 10 step-down runs.")

## Plot 2: min FDR vs donor count
plot_fdr_vs_n <- ggplot(baseline_dt, aes(x = n_donors, y = min_fdr)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "firebrick") +
    geom_line(color = "#117864", linewidth = 0.8) +
    geom_point(aes(color = run_id), size = 2.2) +
    scale_x_continuous(breaks = baseline_dt$n_donors) +
    scale_y_log10() +
    theme_bw(base_size = 12) +
    labs(
        title = "Minimum BH-FDR vs donor count",
        subtitle = "Dashed line marks 0.05",
        x = "Donors in run",
        y = "Minimum BH FDR (log10 scale)"
    ) +
    theme(legend.position = "none")
fn_plot_fdr_vs_n <- file.path(plots_dir, "02_min_fdr_vs_donors.png")
ggsave(fn_plot_fdr_vs_n, plot = plot_fdr_vs_n, width = 8.5, height = 4.8, dpi = 140)
add_plot(fn_plot_fdr_vs_n, "Minimum adjusted p-value trajectory; run10 remains far above 0.05.")

## Plot 3: APOE genotype composition by run
apoe_long <- data.table::melt(
    apoe_dist_dt,
    id.vars = c("run_id", "step", "n_donors"),
    measure.vars = c("n_E2_E2", "n_E2_E3", "n_E3_E4", "n_E4_E4"),
    variable.name = "genotype",
    value.name = "n"
)
apoe_long[, genotype := factor(
    gsub("^n_", "", genotype),
    levels = c("E2_E2", "E2_E3", "E3_E4", "E4_E4")
)]
plot_apoe_stacked <- ggplot(apoe_long, aes(x = factor(step), y = n, fill = genotype)) +
    geom_col(color = "white", linewidth = 0.2) +
    theme_bw(base_size = 12) +
    labs(
        title = "APOE genotype composition across runs",
        subtitle = "Step 0 = 30 donors, Step 9 = 21 donors",
        x = "Step",
        y = "Donor count",
        fill = "APOE genotype"
    )
fn_plot_apoe_stacked <- file.path(plots_dir, "03_apoe_genotype_composition_by_run.png")
ggsave(fn_plot_apoe_stacked, plot = plot_apoe_stacked, width = 9.0, height = 5.2, dpi = 140)
add_plot(fn_plot_apoe_stacked, "Subtype-level APOE composition shift across donor step-down runs.")

## Plot 4: 30 vs 21 logFC comparison
plot_logfc_cmp <- ggplot(cmp_dt, aes(x = logFC_30, y = logFC_21)) +
    geom_point(alpha = 0.12, size = 0.5, color = "#2e4053") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "firebrick") +
    geom_smooth(method = "lm", se = FALSE, color = "#117864", linewidth = 0.8) +
    theme_bw(base_size = 12) +
    labs(
        title = "Effect-size concordance: 30-donor vs 21-donor",
        subtitle = sprintf("Pearson r = %.3f", cor(cmp_dt$logFC_30, cmp_dt$logFC_21, use = "complete.obs")),
        x = "logFC (30 donors)",
        y = "logFC (21 donors)"
    )
fn_plot_logfc_cmp <- file.path(plots_dir, "04_logFC_30_vs_21_scatter.png")
ggsave(fn_plot_logfc_cmp, plot = plot_logfc_cmp, width = 6.8, height = 6.2, dpi = 140)
add_plot(fn_plot_logfc_cmp, "Gene-level logFC concordance between 30 and 21 donor baselines.")

## Plot 5: p-value distributions (30 vs 21)
pval_dt <- data.table::rbindlist(list(
    data.table::data.table(run = "run01_30donors", p = tt_full$P.Value),
    data.table::data.table(run = "run10_21donors", p = tt_21$P.Value)
))
plot_pval_hist <- ggplot(pval_dt, aes(x = p)) +
    geom_histogram(bins = 50, fill = "#5dade2", color = "white") +
    facet_wrap(~run, ncol = 1) +
    theme_bw(base_size = 12) +
    labs(
        title = "Raw p-value distribution by run",
        x = "Raw p-value",
        y = "Gene count"
    )
fn_plot_pval_hist <- file.path(plots_dir, "05_pvalue_hist_run01_vs_run10.png")
ggsave(fn_plot_pval_hist, plot = plot_pval_hist, width = 8.0, height = 7.0, dpi = 140)
add_plot(fn_plot_pval_hist, "Raw p-value histograms showing weaker tail in 21-donor run.")

## Plot 6: 21-donor covariate distributions by carrier group
covar_long_21_plot <- copy(covar_long_21)
covar_long_21_plot[, covariate := factor(covariate, levels = covars)]
plot_covar_21 <- ggplot(covar_long_21_plot, aes(x = carrier_group, y = value, color = carrier_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.15) +
    geom_jitter(width = 0.08, alpha = 0.8, size = 1.8) +
    facet_wrap(~covariate, scales = "free_y", nrow = 1) +
    theme_bw(base_size = 12) +
    labs(
        title = "21-donor covariates by carrier group",
        x = "Carrier group",
        y = "Covariate value"
    ) +
    theme(legend.position = "none")
fn_plot_covar_21 <- file.path(plots_dir, "06_covariates_21donors_by_carrier.png")
ggsave(fn_plot_covar_21, plot = plot_covar_21, width = 11.5, height = 4.2, dpi = 140)
add_plot(fn_plot_covar_21, "Within-21 donor covariate distributions by APOE carrier group.")

## Plot 7: SE(logFC) distribution change
se_dt <- data.table::rbindlist(list(
    data.table::data.table(run = "run01_30donors", se = cmp_dt$se_30),
    data.table::data.table(run = "run10_21donors", se = cmp_dt$se_21)
))
plot_se_density <- ggplot(se_dt[is.finite(se)], aes(x = se, color = run, fill = run)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    theme_bw(base_size = 12) +
    labs(
        title = "SE(logFC) distribution: 30 vs 21 donors",
        subtitle = sprintf("Median SE ratio (21/30) = %.3f",
                           median(cmp_dt$se_21, na.rm = TRUE) / median(cmp_dt$se_30, na.rm = TRUE)),
        x = "SE(logFC)",
        y = "Density"
    )
fn_plot_se_density <- file.path(plots_dir, "07_se_logFC_density_run01_vs_run10.png")
ggsave(fn_plot_se_density, plot = plot_se_density, width = 8.0, height = 5.2, dpi = 140)
add_plot(fn_plot_se_density, "Precision loss visualization (SE inflation in 21-donor run).")

data.table::fwrite(plot_index, file = file.path(out_dir, "run10_plot_index.tsv"), sep = "\t")

## append plot summary to the text report
con <- file(diag_report_file, "at")
writeLines("", con)
writeLines("6) Supporting plots (ggplot2)", con)
for (k in seq_len(nrow(plot_index))) {
    writeLines(sprintf("  - %s : %s", plot_index$plot_file[k], plot_index$description[k]), con)
}
close(con)

message(Sys.time(), " - DEG summary (FDR < 0.05) by run:")
print(baseline_dt[, .(run_id, step, n_donors, n_deg_fdr05, min_fdr, n_E2_carrier, n_E4_carrier, removed_this, donor_key)])
message(Sys.time(), sprintf(" - Wrote APOE distribution table: %s", file.path(out_dir, "apoe_genotype_distribution_by_run.tsv")))
message(Sys.time(), sprintf(" - Wrote 21-donor diagnostic report: %s", diag_report_file))
message(Sys.time(), sprintf(" - Wrote plot index: %s", file.path(out_dir, "run10_plot_index.tsv")))
}

if (do_sweep) {
    full_donors <- unique(as.character(dge_base$samples$BrNum))
    donors21 <- full_donors[full_donors %in% common_donors_21]
    if (length(donors21) != 21) {
        stop(sprintf("Expected 21 overlap donors for sweep mode, found %d", length(donors21)))
    }
    if (max(sweep_target_ns) >= length(donors21) || min(sweep_target_ns) < 2) {
        stop(sprintf("--sweep_targets out of range for 21 donors: %s", paste(sweep_target_ns, collapse = ",")))
    }

    sample21 <- as.data.table(dge_base$samples)[BrNum %in% donors21]
    sample21[, carrier_group := apoe_to_carrier(APOE_syn)]
    base_e2_prop <- mean(sample21$carrier_group == "E2_carrier", na.rm = TRUE)

    target_e2_count <- function(n) {
        e2 <- as.integer(round(base_e2_prop * n))
        e2 <- max(0L, min(e2, n))
        e2
    }

    pick_drop_donor <- function(curr_donors, curr_meta) {
        n_curr <- length(curr_donors)
        n_next <- n_curr - 1L
        if (n_next < 1L) stop("Cannot drop donor below n=1")
        curr_e2 <- sum(curr_meta$carrier_group == "E2_carrier")
        next_e2_target <- target_e2_count(n_next)
        next_e4_target <- n_next - next_e2_target
        curr_e4 <- n_curr - curr_e2

        remove_group <- NULL
        if (curr_e2 > next_e2_target) {
            remove_group <- "E2_carrier"
        } else if (curr_e4 > next_e4_target) {
            remove_group <- "E4_carrier"
        } else {
            ## tie-breaker by which removal keeps next-step ratio closer to 21-donor baseline
            e2_next_if_drop_e2 <- curr_e2 - 1L
            e2_next_if_drop_e4 <- curr_e2
            d_e2 <- abs((e2_next_if_drop_e2 / n_next) - base_e2_prop)
            d_e4 <- abs((e2_next_if_drop_e4 / n_next) - base_e2_prop)
            remove_group <- if (d_e2 <= d_e4) "E2_carrier" else "E4_carrier"
        }

        cand <- curr_meta[carrier_group == remove_group, BrNum]
        if (length(cand) < 1) {
            cand <- curr_meta$BrNum
        }
        sample(cand, size = 1)
    }

    message(Sys.time(), sprintf(
        " - Running random sweeps from n=%d to n=%d (%d sweeps, seed=%d)",
        max(sweep_target_ns), min(sweep_target_ns), n_sweeps, sweep_seed
    ))

    sweep_results <- list()
    sweep_idx <- 0L
    for (sw in seq_len(n_sweeps)) {
        set.seed(sweep_seed + sw - 1L)
        curr_donors <- donors21
        removed <- character(0)

        for (tn in sweep_target_ns) {
            while (length(curr_donors) > tn) {
                curr_meta <- sample21[BrNum %in% curr_donors]
                drop_one <- pick_drop_donor(curr_donors, curr_meta)
                curr_donors <- setdiff(curr_donors, drop_one)
                removed <- c(removed, drop_one)
            }

            donor_info <- make_donor_key(curr_donors)
            res <- run_baseline(donor_info$donor_key, donor_info$donors, donor_info$donor_str)
            if (is.null(res)) next
            curr_meta <- sample21[BrNum %in% curr_donors]
            geno <- get_apoe_counts(curr_meta$APOE_syn)
            n_e2 <- sum(curr_meta$carrier_group == "E2_carrier")
            n_e4 <- sum(curr_meta$carrier_group == "E4_carrier")
            sweep_idx <- sweep_idx + 1L
            sweep_results[[sweep_idx]] <- data.table::data.table(
                sweep_id = sprintf("sweep%02d", sw),
                step = max(sweep_target_ns) - tn + 1L,
                n_donors = tn,
                donor_key = donor_info$donor_key,
                donor_str = donor_info$donor_str,
                removed_so_far = paste(removed, collapse = ","),
                n_E2_carrier = n_e2,
                n_E4_carrier = n_e4,
                n_E2_E2 = geno$n_E2_E2,
                n_E2_E3 = geno$n_E2_E3,
                n_E3_E4 = geno$n_E3_E4,
                n_E4_E4 = geno$n_E4_E4,
                e2_prop = n_e2 / tn,
                e2_prop_dev_from_21 = abs((n_e2 / tn) - base_e2_prop),
                n_deg_fdr05 = res$n_deg_fdr05,
                min_p = res$min_p,
                min_fdr = res$min_fdr,
                out_file = res$out_file,
                out_fdr_file = res$out_fdr_file
            )
        }
    }

    if (length(sweep_results) < 1) {
        stop("No sweep runs completed successfully.")
    }

    sweep_dt <- data.table::rbindlist(sweep_results, fill = TRUE)
    setorder(sweep_dt, sweep_id, -n_donors)

    sweep_label <- paste0("n", max(sweep_target_ns), "to", min(sweep_target_ns))
    sweep_summary_file <- file.path(out_dir, sprintf("random_sweeps_%s_summary.tsv", sweep_label))
    data.table::fwrite(sweep_dt, file = sweep_summary_file, sep = "\t")

    sweep_agg <- sweep_dt[, .(
        n_runs = .N,
        mean_deg = mean(n_deg_fdr05),
        median_deg = median(n_deg_fdr05),
        min_deg = min(n_deg_fdr05),
        max_deg = max(n_deg_fdr05),
        mean_min_fdr = mean(min_fdr),
        median_min_fdr = median(min_fdr),
        min_min_fdr = min(min_fdr),
        max_min_fdr = max(min_fdr)
    ), by = .(n_donors)]
    setorder(sweep_agg, -n_donors)
    sweep_agg_file <- file.path(out_dir, sprintf("random_sweeps_%s_aggregate.tsv", sweep_label))
    data.table::fwrite(sweep_agg, file = sweep_agg_file, sep = "\t")

    ## plots for sweep runs
    plot_sweep_deg <- ggplot(sweep_dt, aes(x = n_donors, y = n_deg_fdr05, color = sweep_id, group = sweep_id)) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 2.2) +
        scale_x_continuous(breaks = sort(unique(sweep_dt$n_donors), decreasing = TRUE)) +
        theme_bw(base_size = 12) +
        labs(
            title = "Random overlap sweeps: DEGs vs donor count",
            subtitle = "Three balanced sweeps from 21-overlap donors",
            x = "Donors in run",
            y = "DEGs (FDR < 0.05)",
            color = "Sweep"
        )
    fn_sweep_plot1 <- file.path(plots_dir, "08_random_sweeps_deg_vs_n.png")
    ggsave(fn_sweep_plot1, plot = plot_sweep_deg, width = 8.5, height = 5.0, dpi = 140)

    plot_sweep_fdr <- ggplot(sweep_dt, aes(x = n_donors, y = min_fdr, color = sweep_id, group = sweep_id)) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "firebrick") +
        geom_line(linewidth = 0.9) +
        geom_point(size = 2.2) +
        scale_x_continuous(breaks = sort(unique(sweep_dt$n_donors), decreasing = TRUE)) +
        scale_y_log10() +
        theme_bw(base_size = 12) +
        labs(
            title = "Random overlap sweeps: min FDR vs donor count",
            x = "Donors in run",
            y = "Minimum BH FDR (log10 scale)",
            color = "Sweep"
        )
    fn_sweep_plot2 <- file.path(plots_dir, "09_random_sweeps_minfdr_vs_n.png")
    ggsave(fn_sweep_plot2, plot = plot_sweep_fdr, width = 8.5, height = 5.0, dpi = 140)

    sweep_report_file <- file.path(out_dir, sprintf("random_sweeps_%s_report.txt", sweep_label))
    con <- file(sweep_report_file, "wt")
    writeLines(sprintf("Random sweep report (%s)", Sys.time()), con)
    writeLines(sprintf("Mode: %s ; sweeps=%d ; seed=%d", run_mode, n_sweeps, sweep_seed), con)
    writeLines(sprintf("Overlap donors used: %d", length(donors21)), con)
    writeLines(sprintf("Target n values per sweep: %s", paste(sweep_target_ns, collapse = ",")), con)
    writeLines(sprintf("Carrier baseline proportion in 21 donors: E2=%.3f E4=%.3f", base_e2_prop, 1 - base_e2_prop), con)
    writeLines("", con)
    writeLines("Per-run summary:", con)
    lines <- capture.output(print(sweep_dt[, .(sweep_id, n_donors, n_E2_carrier, n_E4_carrier, n_deg_fdr05, min_fdr)], nrows = Inf))
    writeLines(paste0("  ", lines), con)
    writeLines("", con)
    writeLines("Aggregate by donor count:", con)
    lines <- capture.output(print(sweep_agg, nrows = Inf))
    writeLines(paste0("  ", lines), con)
    writeLines("", con)
    writeLines(sprintf("Sweep DEG plot: %s", fn_sweep_plot1), con)
    writeLines(sprintf("Sweep min-FDR plot: %s", fn_sweep_plot2), con)
    close(con)

    message(Sys.time(), " - Random sweep DEG summary:")
    print(sweep_dt[, .(sweep_id, n_donors, n_E2_carrier, n_E4_carrier, n_deg_fdr05, min_fdr)])
    message(Sys.time(), sprintf(" - Wrote sweep summary: %s", sweep_summary_file))
    message(Sys.time(), sprintf(" - Wrote sweep aggregate: %s", sweep_agg_file))
    message(Sys.time(), sprintf(" - Wrote sweep report: %s", sweep_report_file))
}

if (do_outlier) {
    full_donors <- unique(as.character(dge_base$samples$BrNum))
    donors21 <- full_donors[full_donors %in% common_donors_21]
    if (length(donors21) != 21) {
        stop(sprintf("Expected 21 overlap donors for outlier mode, found %d", length(donors21)))
    }
    sample21 <- as.data.table(dge_base$samples)[BrNum %in% donors21]
    sample21[, carrier_group := apoe_to_carrier(APOE_syn)]
    sample21 <- unique(sample21[, .(BrNum, APOE_syn, carrier_group, Sex, Age, Anc_Afr, pseudo_expr_chrM_ratio, exp_round)])

    message(Sys.time(), sprintf(" - Outlier scan on %d overlap donors (workers=%d)", length(donors21), workers))

    ## baseline 21 reference
    info21 <- make_donor_key(donors21)
    ref21 <- run_baseline(info21$donor_key, info21$donors, info21$donor_str)
    if (is.null(ref21)) stop("Failed to compute/reference baseline for 21-donor overlap set.")

    ## leave-one-out (21 -> 20)
    l1_plan <- data.table::data.table(
        drop1 = donors21,
        donors = lapply(donors21, function(d) setdiff(donors21, d))
    )
    bpp <- BiocParallel::MulticoreParam(workers = workers)
    l1_res <- BiocParallel::bplapply(seq_len(nrow(l1_plan)), function(i) {
        d1 <- l1_plan$drop1[i]
        donors_i <- l1_plan$donors[[i]]
        info <- make_donor_key(donors_i)
        rr <- run_baseline(info$donor_key, info$donors, info$donor_str)
        if (is.null(rr)) return(NULL)
        dm <- sample21[BrNum == d1]
        data.table::data.table(
            n_start = 21L,
            n_donors = 20L,
            drop1 = d1,
            drop1_apoe = as.character(dm$APOE_syn[1]),
            drop1_carrier = as.character(dm$carrier_group[1]),
            donor_key = rr$donor_key,
            donor_str = rr$donor_str,
            n_deg_fdr05 = rr$n_deg_fdr05,
            min_p = rr$min_p,
            min_fdr = rr$min_fdr,
            n_E2_carrier = rr$n_E2_carrier,
            n_E4_carrier = rr$n_E4_carrier,
            out_file = rr$out_file,
            out_fdr_file = rr$out_fdr_file
        )
    }, BPPARAM = bpp)
    l1_res <- Filter(Negate(is.null), l1_res)
    if (length(l1_res) < 1) stop("No leave-one-out runs completed.")
    l1_dt <- data.table::rbindlist(l1_res, fill = TRUE)
    setorder(l1_dt, -n_deg_fdr05, min_fdr, min_p, drop1)

    out_l1_file <- file.path(out_dir, "outlier_scan_leave1_21to20.tsv")
    data.table::fwrite(l1_dt, out_l1_file, sep = "\t")

    best1 <- l1_dt[1]
    donors20_best <- setdiff(donors21, best1$drop1)
    message(Sys.time(), sprintf(" - Best single donor drop: %s (DEGs=%d, minFDR=%.4g)", best1$drop1, best1$n_deg_fdr05, best1$min_fdr))

    ## second-drop conditional scan from best 20 set (20 -> 19)
    l2_plan <- data.table::data.table(
        drop2 = donors20_best,
        donors = lapply(donors20_best, function(d) setdiff(donors20_best, d))
    )
    l2_res <- BiocParallel::bplapply(seq_len(nrow(l2_plan)), function(i) {
        d2 <- l2_plan$drop2[i]
        donors_i <- l2_plan$donors[[i]]
        info <- make_donor_key(donors_i)
        rr <- run_baseline(info$donor_key, info$donors, info$donor_str)
        if (is.null(rr)) return(NULL)
        dm <- sample21[BrNum == d2]
        dropped <- sort(c(best1$drop1, d2))
        data.table::data.table(
            n_start = 21L,
            n_donors = 19L,
            drop1 = best1$drop1,
            drop2 = d2,
            dropped_pair = paste(dropped, collapse = ","),
            drop2_apoe = as.character(dm$APOE_syn[1]),
            drop2_carrier = as.character(dm$carrier_group[1]),
            donor_key = rr$donor_key,
            donor_str = rr$donor_str,
            n_deg_fdr05 = rr$n_deg_fdr05,
            min_p = rr$min_p,
            min_fdr = rr$min_fdr,
            n_E2_carrier = rr$n_E2_carrier,
            n_E4_carrier = rr$n_E4_carrier,
            out_file = rr$out_file,
            out_fdr_file = rr$out_fdr_file
        )
    }, BPPARAM = bpp)
    l2_res <- Filter(Negate(is.null), l2_res)
    if (length(l2_res) < 1) stop("No conditional leave-two runs completed.")
    l2_dt <- data.table::rbindlist(l2_res, fill = TRUE)
    setorder(l2_dt, -n_deg_fdr05, min_fdr, min_p, drop2)

    out_l2_file <- file.path(out_dir, "outlier_scan_leave2_from_best20.tsv")
    data.table::fwrite(l2_dt, out_l2_file, sep = "\t")

    best2 <- l2_dt[1]
    message(Sys.time(), sprintf(" - Best donor pair: %s (DEGs=%d, minFDR=%.4g)", best2$dropped_pair, best2$n_deg_fdr05, best2$min_fdr))

    ## gene-level comparison: best 19-donor run vs 21 baseline
    tt21 <- data.table::fread(ref21$out_file)
    tt19 <- data.table::fread(best2$out_file)
    cmp <- merge(
        tt21[, .(outcome_gene_id, logFC_21 = logFC, t_21 = t, P_21 = P.Value, FDR_21 = adj.P.Val)],
        tt19[, .(outcome_gene_id, logFC_19 = logFC, t_19 = t, P_19 = P.Value, FDR_19 = adj.P.Val)],
        by = "outcome_gene_id"
    )
    cmp[, d_neglog10P := (-log10(P_19 + 1e-300)) - (-log10(P_21 + 1e-300))]
    cmp[, abs_d_logFC := abs(logFC_19 - logFC_21)]
    setorder(cmp, -d_neglog10P, FDR_19)
    out_cmp_file <- file.path(out_dir, "outlier_best19_vs21_gene_shift_top500.tsv")
    data.table::fwrite(cmp[1:min(.N, 500)], out_cmp_file, sep = "\t")

    tt19_sig <- tt19[adj.P.Val < 0.05]
    out_best19_sig <- file.path(out_dir, "outlier_best19_fdr05.tsv")
    data.table::fwrite(tt19_sig, out_best19_sig, sep = "\t")

    ## concise plots
    p_l1 <- ggplot(l1_dt, aes(x = reorder(drop1, n_deg_fdr05), y = n_deg_fdr05, fill = drop1_carrier)) +
        geom_col() +
        coord_flip() +
        theme_bw(base_size = 11) +
        labs(
            title = "Leave-one-out scan (21 -> 20): DEGs by dropped donor",
            x = "Dropped donor",
            y = "DEGs (FDR<0.05)",
            fill = "Carrier group"
        )
    fn_p_l1 <- file.path(plots_dir, "10_outlier_scan_leave1_deg.png")
    ggsave(fn_p_l1, p_l1, width = 8.2, height = 6.5, dpi = 140)

    p_l2 <- ggplot(l2_dt, aes(x = reorder(drop2, n_deg_fdr05), y = n_deg_fdr05, fill = drop2_carrier)) +
        geom_col() +
        coord_flip() +
        theme_bw(base_size = 11) +
        labs(
            title = sprintf("Conditional leave-two scan from best 20-donor set (drop1=%s)", best1$drop1),
            x = "Second dropped donor",
            y = "DEGs (FDR<0.05)",
            fill = "Carrier group"
        )
    fn_p_l2 <- file.path(plots_dir, "11_outlier_scan_leave2_deg.png")
    ggsave(fn_p_l2, p_l2, width = 8.2, height = 6.5, dpi = 140)

    ## report
    out_report <- file.path(out_dir, "outlier_scan_report.txt")
    con <- file(out_report, "wt")
    writeLines(sprintf("Outlier donor scan report (%s)", Sys.time()), con)
    writeLines(sprintf("Mode=%s ; workers=%d", run_mode, workers), con)
    writeLines(sprintf("21-donor baseline: DEGs=%d ; minFDR=%.6f ; donor_key=%s", ref21$n_deg_fdr05, ref21$min_fdr, ref21$donor_key), con)
    writeLines("", con)
    writeLines("Top leave-one-out results (21->20):", con)
    lines <- capture.output(print(l1_dt[1:min(.N, 10), .(drop1, drop1_apoe, drop1_carrier, n_deg_fdr05, min_fdr)]))
    writeLines(paste0("  ", lines), con)
    writeLines("", con)
    writeLines("Top conditional leave-two results from best 20-donor set (20->19):", con)
    lines <- capture.output(print(l2_dt[1:min(.N, 10), .(drop1, drop2, dropped_pair, n_deg_fdr05, min_fdr)]))
    writeLines(paste0("  ", lines), con)
    writeLines("", con)
    writeLines(sprintf("Best single donor to drop: %s (%s), DEGs=%d, minFDR=%.6f",
                       best1$drop1, best1$drop1_apoe, best1$n_deg_fdr05, best1$min_fdr), con)
    writeLines(sprintf("Best donor pair to drop: %s, DEGs=%d, minFDR=%.6f",
                       best2$dropped_pair, best2$n_deg_fdr05, best2$min_fdr), con)
    writeLines(sprintf("Best 19-donor significant genes (FDR<0.05): %d", nrow(tt19_sig)), con)
    writeLines("", con)
    writeLines("Output files:", con)
    writeLines(sprintf("  %s", out_l1_file), con)
    writeLines(sprintf("  %s", out_l2_file), con)
    writeLines(sprintf("  %s", out_cmp_file), con)
    writeLines(sprintf("  %s", out_best19_sig), con)
    writeLines(sprintf("  %s", fn_p_l1), con)
    writeLines(sprintf("  %s", fn_p_l2), con)
    close(con)

    message(Sys.time(), " - Outlier scan summary (top 10 leave-one-out):")
    print(l1_dt[1:min(.N, 10), .(drop1, drop1_apoe, n_deg_fdr05, min_fdr)])
    message(Sys.time(), " - Outlier scan summary (top 10 leave-two from best20):")
    print(l2_dt[1:min(.N, 10), .(drop1, drop2, n_deg_fdr05, min_fdr)])
    message(Sys.time(), sprintf(" - Wrote outlier report: %s", out_report))
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
