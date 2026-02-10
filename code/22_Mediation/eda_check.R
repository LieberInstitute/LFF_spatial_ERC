#!/bin/env Rscript

## Focused EDA for donor-level impact on 21->20 overlap baseline DGE
## Default focus donor: Br5529
##
## Purpose:
##   Quantify *why* dropping one donor (default Br5529) restores DEG discovery.
##   This script is intentionally verbose and auditable: it logs each step,
##   writes intermediate/summary tables, and produces diagnostics plots.
##
## Key idea:
##   The question is not only "does leave-one-out change DEGs?" but whether the
##   focus donor is atypical in covariates and/or expression space, and whether
##   its expression profile specifically attenuates the carrier contrast.

library("data.table")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("SummarizedExperiment")
library("here")
library("getopt")
library("ggplot2")
library("msigdbr")
library("fgsea")

## Map 4-level APOE genotype into carrier groups used for the main contrast.
apoe_to_carrier <- function(x) {
    x <- as.character(x)
    out <- ifelse(x %in% c("E2.E2", "E2.E3"), "E2_carrier",
                  ifelse(x %in% c("E3.E4", "E4.E4"), "E4_carrier", NA_character_))
    factor(out, levels = c("E2_carrier", "E4_carrier"))
}

## Convenience hash helper for deterministic run identifiers if needed later.
md5_string <- function(x) {
    tf <- tempfile()
    writeLines(x, con = tf, useBytes = TRUE)
    md5 <- as.character(tools::md5sum(tf))
    unlink(tf)
    md5
}

## Consistent numeric formatting for report text.
fmt <- function(x, digits = 3) {
    if (length(x) < 1 || is.na(x) || !is.finite(x)) return("NA")
    formatC(x, format = "f", digits = digits)
}

## CLI options
##   - focus_donor: donor to investigate
##   - out_dir: where all tables/plots/report/log are written
##   - top_n: number of top genes to keep in detailed tables
##   - cohort_label: cohort stamp embedded in plot/report titles
##   - plot_profile: full (all plots) or explanatory (08-11 only)
##   - assert_n_overlap: hard-check overlap donor count before analysis
##   - highlight_donors: comma-separated donors to always label in donor-granular plots
opt_spec <- matrix(
    c(
        "focus_donor", "d", "1", "character", "Donor to investigate (default=Br5529)",
        "out_dir", "o", "1", "character", "Output directory (default=code/22_Mediation/out-dge-test/eda_br5529)",
        "top_n", "n", "1", "integer", "Top N genes to include in detailed tables (default=200)",
        "cohort_label", "c", "1", "character", "Label for cohort stamp (default=21-overlap)",
        "plot_profile", "p", "1", "character", "Plot profile: full|explanatory (default=full)",
        "assert_n_overlap", "a", "1", "integer", "Assert expected overlap donor count (default=21)",
        "highlight_donors", "h", "1", "character", "Comma-separated donor IDs to always label (default=Br6423,Br6263)"
    ),
    ncol = 5, byrow = TRUE
)
opt <- tryCatch(getopt(opt_spec), error = function(e) list())
focus_donor <- if (is.null(opt$focus_donor) || is.na(opt$focus_donor)) "Br5529" else as.character(opt$focus_donor)
out_dir <- if (is.null(opt$out_dir) || is.na(opt$out_dir)) {
    here("code", "22_Mediation", "out-dge-test", "eda_br5529")
} else {
    as.character(opt$out_dir)
}
top_n <- if (is.null(opt$top_n) || is.na(opt$top_n)) 200L else as.integer(opt$top_n)
if (is.na(top_n) || top_n < 10L) top_n <- 200L
cohort_label <- if (is.null(opt$cohort_label) || is.na(opt$cohort_label)) "21-overlap" else as.character(opt$cohort_label)
plot_profile <- if (is.null(opt$plot_profile) || is.na(opt$plot_profile)) "full" else tolower(as.character(opt$plot_profile))
assert_n_overlap <- if (is.null(opt$assert_n_overlap) || is.na(opt$assert_n_overlap)) 21L else as.integer(opt$assert_n_overlap)
highlight_donors_str <- if (is.null(opt$highlight_donors) || is.na(opt$highlight_donors)) "Br6423,Br6263" else as.character(opt$highlight_donors)
highlight_donors <- unique(trimws(strsplit(highlight_donors_str, ",", fixed = TRUE)[[1]]))
highlight_donors <- highlight_donors[nchar(highlight_donors) > 0]
if (!(plot_profile %in% c("full", "explanatory"))) stop("Invalid --plot_profile. Use full or explanatory.")
if (is.na(assert_n_overlap) || assert_n_overlap < 2L) stop("Invalid --assert_n_overlap, expected >=2.")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, sprintf("eda_check_%s.log", focus_donor))
if (file.exists(log_file)) invisible(file.remove(log_file))
log_msg <- function(...) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", paste0(..., collapse = ""))
    cat(line, "\n")
    cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Starting EDA check for focus donor: ", focus_donor)

## Paths
pb_fn <- here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", "sce_pseudo_DGE-cell_type_anno.RDS")
f_lca_pbcounts <- here("processed-data", "22_Mediation", "lc_astro_pseudobulk_counts_by_donor.rds")

## Load data
##   - ERC Oligo.3 pseudobulk input
##   - LC pseudobulk donor list used only to define the 21-donor overlap set
log_msg("Loading ERC pseudobulk and LC donor overlap data")
sce_pb <- readRDS(pb_fn)
sce_pb <- sce_pb[, sce_pb$BrNum != "Br1289"]
sce_base <- sce_pb[, sce_pb$registration_variable == "Oligo.3"]

lca_pb <- readRDS(f_lca_pbcounts)
lca_brnums <- colnames(lca_pb)

meta_all <- as.data.table(as.data.frame(colData(sce_base)))
meta_all <- unique(meta_all[, .(BrNum, APOE_syn, Sex, Age, Anc_Afr, pseudo_expr_chrM_ratio, exp_round)])
overlap_donors <- sort(intersect(lca_brnums, as.character(meta_all$BrNum)))
overlap_donor_hash <- md5_string(paste(overlap_donors, collapse = ","))

if (!(focus_donor %in% overlap_donors)) {
    stop(sprintf("Focus donor %s is not in the 21 overlap donors.", focus_donor))
}
if (length(overlap_donors) != assert_n_overlap) {
    stop(sprintf("Overlap donor assertion failed: observed=%d expected=%d", length(overlap_donors), assert_n_overlap))
}

log_msg("Overlap donors: ", length(overlap_donors))
log_msg("Focus donor in overlap set: ", focus_donor)
log_msg("Overlap donor hash: ", overlap_donor_hash)
log_msg("Plot profile: ", plot_profile, " ; cohort label: ", cohort_label)
log_msg("Always-highlight donors: ", paste(highlight_donors, collapse = ","))

## Model helper (same baseline logic as baseline_testing.R)
## This ensures EDA is aligned with the exact DE workflow used in screening.
baseline_formula <- as.formula("~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio")
carrier_contrast_expr <- "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)"

run_model <- function(donors) {
    donors <- sort(unique(as.character(donors)))
    out_brnum <- as.character(sce_base$BrNum)
    idx <- match(donors, out_brnum)
    idx <- idx[!is.na(idx)]
    sce_sub <- sce_base[, idx]

    dge <- edgeR::calcNormFactors(sce_sub)
    sample_df <- as.data.frame(dge$samples)
    sample_df$APOE_syn <- droplevels(factor(sample_df$APOE_syn))
    sample_df$Sex <- droplevels(factor(sample_df$Sex))
    sample_df$exp_round <- droplevels(factor(sample_df$exp_round))

    des <- model.matrix(baseline_formula, data = sample_df)
    colnames(des) <- gsub(colnames(des), pattern = "_syn", replacement = "_")

    ## Recompute filtering and normalization for each donor set independently.
    keep <- edgeR::filterByExpr(dge, design = des)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)

    ## Fit with block + sample-weights, matching production pipeline behavior.
    fit <- voomLmFit(
        dge,
        design = des,
        block = as.factor(sample_df$exp_round),
        adaptive.span = TRUE,
        sample.weights = TRUE
    )

    cont <- makeContrasts(
        contrasts = carrier_contrast_expr,
        levels = des
    )
    colnames(cont) <- "carrier"
    fit2 <- eBayes(contrasts.fit(fit, cont))
    tt <- as.data.table(topTable(fit2, coef = 1, number = Inf, adjust.method = "BH"), keep.rownames = "gene")

    list(
        donors = donors,
        dge = dge,
        sample_df = sample_df,
        design = des,
        fit = fit,
        tt = tt
    )
}

## Fit the paired models:
##   1) 21-overlap baseline
##   2) 20-donor model after dropping focus donor
## This directly quantifies the effect size of donor exclusion.
log_msg("Running 21-donor model")
res21 <- run_model(overlap_donors)
donors20 <- setdiff(overlap_donors, focus_donor)
if (length(donors20) != (length(overlap_donors) - 1L) || !setequal(donors20, setdiff(overlap_donors, focus_donor))) {
    stop("Cohort integrity assertion failed for 20-donor set.")
}
log_msg("Running 20-donor model (drop ", focus_donor, ")")
res20 <- run_model(donors20)

n_deg21 <- sum(res21$tt$adj.P.Val < 0.05, na.rm = TRUE)
n_deg20 <- sum(res20$tt$adj.P.Val < 0.05, na.rm = TRUE)
min_fdr21 <- min(res21$tt$adj.P.Val, na.rm = TRUE)
min_fdr20 <- min(res20$tt$adj.P.Val, na.rm = TRUE)

log_msg("Results: 21-donor DEGs=", n_deg21, " (minFDR=", signif(min_fdr21, 4), "); 20-donor DEGs=", n_deg20, " (minFDR=", signif(min_fdr20, 4), ")")

## Donor-level diagnostics (21-donor run)
## All outlier metrics below are computed in the 21-donor baseline space.
diag_dt <- copy(as.data.table(res21$sample_df)[, .(
    BrNum = as.character(BrNum),
    APOE_syn = as.character(APOE_syn),
    Sex = as.character(Sex),
    Age = as.numeric(Age),
    Anc_Afr = as.numeric(Anc_Afr),
    pseudo_expr_chrM_ratio = as.numeric(pseudo_expr_chrM_ratio),
    exp_round = as.character(exp_round)
)])
diag_dt <- unique(diag_dt, by = "BrNum")
diag_dt[, carrier := as.character(apoe_to_carrier(APOE_syn))]

## Covariate z-scores:
##   - global z across 21 donors
##   - within-carrier z (helps separate carrier composition from donor-specific shift)
num_cov <- c("Age", "Anc_Afr", "pseudo_expr_chrM_ratio")
for (v in num_cov) {
    z <- as.numeric(scale(diag_dt[[v]]))
    diag_dt[[paste0(v, "_z")]] <- z
}
for (v in num_cov) {
    diag_dt[[paste0(v, "_z_within_carrier")]] <- diag_dt[, (get(v) - mean(get(v), na.rm = TRUE)) / sd(get(v), na.rm = TRUE), by = carrier]$V1
}

## Mahalanobis distance on numeric model covariates.
## Interpretation: high values indicate multivariate covariate atypicality.
cov_mat <- as.matrix(scale(diag_dt[, ..num_cov]))
data.table::set(diag_dt, j = "mahalanobis_cov", value = mahalanobis(cov_mat, center = colMeans(cov_mat), cov = cov(cov_mat)))

## Design leverage h_ii (from X'X).
## Interpretation: high leverage means donor has strong influence on fitted coefficients.
X <- as.matrix(res21$design)
XtX_inv <- solve(crossprod(X))
hii <- rowSums((X %*% XtX_inv) * X)
lev_dt <- data.table(BrNum = as.character(res21$sample_df$BrNum), leverage = hii)
diag_dt <- merge(diag_dt, lev_dt, by = "BrNum", all.x = TRUE)

## Voom-level quality diagnostics:
##   - mean observational weight (lower may suggest noisier sample)
##   - approximate residual MAD after linear fit (higher residual burden = less model-conforming)
obs_w <- res21$fit$EList$weights
if (is.null(obs_w) || !is.matrix(obs_w)) {
    stop("Expected observational weights matrix in fit$EList$weights.")
}
sample_weight_mean <- colMeans(obs_w, na.rm = TRUE)
diag_dt <- merge(
    diag_dt,
    data.table(BrNum = as.character(res21$sample_df$BrNum), mean_obs_weight = as.numeric(sample_weight_mean)),
    by = "BrNum", all.x = TRUE
)

E <- res21$fit$EList$E
B <- res21$fit$coefficients
fitted_E <- B %*% t(X)
resid_E <- E - fitted_E
sample_resid_mad <- apply(abs(resid_E), 2, median, na.rm = TRUE)
diag_dt <- merge(
    diag_dt,
    data.table(BrNum = as.character(res21$sample_df$BrNum), approx_resid_mad = as.numeric(sample_resid_mad)),
    by = "BrNum", all.x = TRUE
)

## Expression-space diagnostics:
##   - PCA distance from cohort center (PC1-3)
##   - median sample-sample Spearman correlation
##   - burden of extreme genes (|z| > 3 across donors)
logcpm21 <- edgeR::cpm(res21$dge, log = TRUE, prior.count = 2)
colnames(logcpm21) <- as.character(res21$sample_df$BrNum)

pc <- prcomp(t(logcpm21), center = TRUE, scale. = FALSE)
pc123_center <- colMeans(pc$x[, 1:3, drop = FALSE])
pc_dist <- sqrt(rowSums((pc$x[, 1:3, drop = FALSE] - matrix(pc123_center, nrow = nrow(pc$x), ncol = 3, byrow = TRUE))^2))
diag_dt <- merge(diag_dt, data.table(BrNum = rownames(pc$x), pca_dist_pc123 = as.numeric(pc_dist)), by = "BrNum", all.x = TRUE)

cor_spearman <- cor(logcpm21, method = "spearman")
med_cor <- sapply(colnames(cor_spearman), function(s) median(cor_spearman[s, setdiff(colnames(cor_spearman), s)], na.rm = TRUE))
diag_dt <- merge(diag_dt, data.table(BrNum = names(med_cor), median_spearman_corr = as.numeric(med_cor)), by = "BrNum", all.x = TRUE)

z_gene <- t(scale(t(logcpm21)))
n_extreme <- colSums(abs(z_gene) > 3, na.rm = TRUE)
diag_dt <- merge(diag_dt, data.table(BrNum = names(n_extreme), n_extreme_genes_absz_gt3 = as.integer(n_extreme)), by = "BrNum", all.x = TRUE)

## Rank diagnostics.
## All rank_* fields use rank=1 as "most outlying in that direction".
rank_high <- function(x) rank(-x, ties.method = "min", na.last = "keep")
rank_low <- function(x) rank(x, ties.method = "min", na.last = "keep")
diag_dt[, mahalanobis_rank_high := rank_high(mahalanobis_cov)]
diag_dt[, leverage_rank_high := rank_high(leverage)]
diag_dt[, pca_dist_rank_high := rank_high(pca_dist_pc123)]
diag_dt[, residual_burden_rank_high := rank_high(approx_resid_mad)]
diag_dt[, low_corr_rank := rank_low(median_spearman_corr)]
diag_dt[, extreme_gene_burden_rank_high := rank_high(n_extreme_genes_absz_gt3)]
diag_dt[, chrM_ratio_rank_high := rank_high(pseudo_expr_chrM_ratio)]

## Donors to label on donor-granular plots:
##  - focus donor,
##  - user-requested always-highlight donors,
##  - "outlier-like" donors by top ranks on key metrics.
diag_dt[, outlier_like := (pca_dist_rank_high <= 3L) | (low_corr_rank <= 3L) | (extreme_gene_burden_rank_high <= 3L) | (mahalanobis_rank_high <= 3L)]
diag_dt[, near_outlier_like := (pca_dist_rank_high <= 5L) | (low_corr_rank <= 5L) | (extreme_gene_burden_rank_high <= 5L) | (mahalanobis_rank_high <= 5L)]
auto_label_donors <- unique(c(
    diag_dt[outlier_like == TRUE, BrNum],
    diag_dt[near_outlier_like == TRUE, BrNum]
))
label_donors <- sort(unique(c(focus_donor, highlight_donors, auto_label_donors)))
label_donors <- label_donors[label_donors %in% diag_dt$BrNum]
diag_dt[, label_group := fifelse(
    BrNum == focus_donor,
    "focus",
    fifelse(
        BrNum %in% highlight_donors,
        "requested",
        fifelse(outlier_like == TRUE, "outlier", fifelse(near_outlier_like == TRUE, "near-outlier", "other"))
    )
)]
diag_dt[, label_group := factor(label_group, levels = c("focus", "requested", "outlier", "near-outlier", "other"))]

## Donor labels for axes:
##  - "<- focus": active donor under test
##  - "[req]": user-requested donor labels
##  - "[out]": top outlier-like donors
##  - "[near]": near-outlier donors
donor_axis_labels <- setNames(as.character(diag_dt$BrNum), as.character(diag_dt$BrNum))
for (d in names(donor_axis_labels)) {
    g <- as.character(diag_dt[BrNum == d, label_group][1])
    suffix <- if (g == "focus") {
        " <- focus"
    } else if (g == "requested") {
        " [req]"
    } else if (g == "outlier") {
        " [out]"
    } else if (g == "near-outlier") {
        " [near]"
    } else {
        ""
    }
    donor_axis_labels[d] <- paste0(donor_axis_labels[d], suffix)
}

## Gene-level attenuation diagnostics for focus donor.
## We evaluate whether focus donor tends to be shifted within its own carrier group
## and whether removing it systematically increases the carrier contrast.
focus_carrier <- diag_dt[BrNum == focus_donor, carrier][1]
group_ids <- diag_dt[carrier == focus_carrier, BrNum]
if (!(focus_donor %in% group_ids)) {
    stop("Focus donor not present in its own carrier group.")
}
group_ids_wo <- setdiff(group_ids, focus_donor)
opp_ids <- diag_dt[carrier != focus_carrier, BrNum]

focus_expr <- logcpm21[, focus_donor]
group_mean <- rowMeans(logcpm21[, group_ids, drop = FALSE], na.rm = TRUE)
group_sd <- apply(logcpm21[, group_ids, drop = FALSE], 1, sd, na.rm = TRUE)
focus_z_in_group <- (focus_expr - group_mean) / group_sd

group_mean_wo <- rowMeans(logcpm21[, group_ids_wo, drop = FALSE], na.rm = TRUE)
opp_mean <- rowMeans(logcpm21[, opp_ids, drop = FALSE], na.rm = TRUE)

contrast_with_focus <- opp_mean - group_mean
contrast_without_focus <- opp_mean - group_mean_wo
contrast_gain_drop_focus <- contrast_without_focus - contrast_with_focus

tt21 <- res21$tt[, .(gene, logFC_21 = logFC, P_21 = P.Value, FDR_21 = adj.P.Val)]
tt20 <- res20$tt[, .(gene, logFC_20 = logFC, P_20 = P.Value, FDR_20 = adj.P.Val)]

gene_dt <- data.table(
    gene = rownames(logcpm21),
    focus_expr = as.numeric(focus_expr),
    focus_z_in_group = as.numeric(focus_z_in_group),
    contrast_with_focus = as.numeric(contrast_with_focus),
    contrast_without_focus = as.numeric(contrast_without_focus),
    contrast_gain_drop_focus = as.numeric(contrast_gain_drop_focus)
)
gene_dt <- merge(gene_dt, tt21, by = "gene", all.x = TRUE)
gene_dt <- merge(gene_dt, tt20, by = "gene", all.x = TRUE)
gene_dt[, is_deg20 := FDR_20 < 0.05]
gene_dt[, is_deg21 := FDR_21 < 0.05]

## Summary stats on focus z pattern
z_all <- gene_dt$focus_z_in_group[is.finite(gene_dt$focus_z_in_group)]
z_deg20 <- gene_dt[is_deg20 == TRUE & is.finite(focus_z_in_group), focus_z_in_group]
gain_deg20 <- gene_dt[is_deg20 == TRUE & is.finite(contrast_gain_drop_focus), contrast_gain_drop_focus]

## One-sided sign tests: is the fraction of positive z larger than expected by chance?
sign_test_all <- binom.test(sum(z_all > 0), length(z_all), p = 0.5, alternative = "greater")
sign_test_deg20 <- if (length(z_deg20) > 0) binom.test(sum(z_deg20 > 0), length(z_deg20), p = 0.5, alternative = "greater") else NULL

summary_dt <- data.table(
    metric = c(
        "n_genes_tested_21",
        "n_deg_21_fdr05",
        "n_deg_20_fdr05_drop_focus",
        "min_fdr_21",
        "min_fdr_20_drop_focus",
        "focus_z_all_median",
        "focus_z_all_mean",
        "focus_z_all_frac_gt0",
        "focus_z_all_frac_gt1.5",
        "focus_z_all_frac_lt_minus1.5",
        "focus_z_deg20_median",
        "focus_z_deg20_mean",
        "focus_z_deg20_frac_gt0",
        "focus_z_deg20_frac_gt1.5",
        "contrast_gain_deg20_median",
        "contrast_gain_deg20_mean",
        "contrast_gain_deg20_frac_gt0",
        "binom_p_all_z_gt0",
        "binom_p_deg20_z_gt0"
    ),
    value = c(
        nrow(gene_dt),
        n_deg21,
        n_deg20,
        min_fdr21,
        min_fdr20,
        median(z_all, na.rm = TRUE),
        mean(z_all, na.rm = TRUE),
        mean(z_all > 0, na.rm = TRUE),
        mean(z_all > 1.5, na.rm = TRUE),
        mean(z_all < -1.5, na.rm = TRUE),
        if (length(z_deg20) > 0) median(z_deg20, na.rm = TRUE) else NA_real_,
        if (length(z_deg20) > 0) mean(z_deg20, na.rm = TRUE) else NA_real_,
        if (length(z_deg20) > 0) mean(z_deg20 > 0, na.rm = TRUE) else NA_real_,
        if (length(z_deg20) > 0) mean(z_deg20 > 1.5, na.rm = TRUE) else NA_real_,
        if (length(gain_deg20) > 0) median(gain_deg20, na.rm = TRUE) else NA_real_,
        if (length(gain_deg20) > 0) mean(gain_deg20, na.rm = TRUE) else NA_real_,
        if (length(gain_deg20) > 0) mean(gain_deg20 > 0, na.rm = TRUE) else NA_real_,
        as.numeric(sign_test_all$p.value),
        if (!is.null(sign_test_deg20)) as.numeric(sign_test_deg20$p.value) else NA_real_
    )
)

## Top gene tables for audit:
##   - genes that are DE in 20-donor model
##   - genes with largest carrier-contrast gain after dropping focus donor
top_deg20 <- gene_dt[is_deg20 == TRUE][order(FDR_20, -abs(focus_z_in_group))][1:min(.N, top_n)]
top_by_gain <- gene_dt[is.finite(contrast_gain_drop_focus)][order(-contrast_gain_drop_focus)][1:min(.N, top_n)]

## Within-focus-carrier comparison on the same DEG set.
## If focus donor ranks highest here, that supports donor-specific attenuation.
carrier_cmp_dt <- data.table()
if (sum(gene_dt$is_deg20, na.rm = TRUE) > 0) {
    deg20_genes <- gene_dt[is_deg20 == TRUE, gene]
    zmat_carrier <- t(scale(t(logcpm21[deg20_genes, group_ids, drop = FALSE])))
    carrier_cmp_dt <- data.table(
        BrNum = colnames(zmat_carrier),
        median_z_deg20 = apply(zmat_carrier, 2, median, na.rm = TRUE),
        mean_z_deg20 = apply(zmat_carrier, 2, mean, na.rm = TRUE),
        frac_z_gt0_deg20 = apply(zmat_carrier, 2, function(x) mean(x > 0, na.rm = TRUE)),
        frac_z_gt1p5_deg20 = apply(zmat_carrier, 2, function(x) mean(x > 1.5, na.rm = TRUE))
    )
    setorder(carrier_cmp_dt, -median_z_deg20)
    carrier_cmp_dt[, median_z_rank_high := rank(-median_z_deg20, ties.method = "min")]
}

## Write tables
f_diag <- file.path(out_dir, "donor_diagnostics_21donors.tsv")
f_focus <- file.path(out_dir, sprintf("focus_donor_%s_diagnostics.tsv", focus_donor))
f_gene <- file.path(out_dir, "gene_level_focus_effect.tsv.gz")
f_sum <- file.path(out_dir, "focus_effect_summary.tsv")
f_top_deg20 <- file.path(out_dir, sprintf("top_deg20_focus_%s.tsv", focus_donor))
f_top_gain <- file.path(out_dir, sprintf("top_genes_by_contrast_gain_drop_%s.tsv", focus_donor))
f_carrier_cmp <- file.path(out_dir, sprintf("carrier_group_deg20_z_comparison_%s.tsv", focus_donor))
f_l1_influence <- file.path(out_dir, sprintf("leave1_influence_21to20_%s.tsv", focus_donor))
f_fgsea <- file.path(out_dir, sprintf("high_gain_fgsea_%s.tsv", focus_donor))
f_plot_manifest <- file.path(out_dir, "plot_manifest.tsv")
f_scope <- file.path(out_dir, "readme_scope.txt")
f_digest <- file.path(out_dir, sprintf("eda_report_%s_digest.md", focus_donor))
f_slide_md <- file.path(out_dir, sprintf("%s_one_slide.md", tolower(focus_donor)))

fwrite(diag_dt[order(low_corr_rank, -pca_dist_pc123)], f_diag, sep = "\t")
fwrite(diag_dt[BrNum == focus_donor], f_focus, sep = "\t")
fwrite(gene_dt, f_gene, sep = "\t")
fwrite(summary_dt, f_sum, sep = "\t")
fwrite(top_deg20, f_top_deg20, sep = "\t")
fwrite(top_by_gain, f_top_gain, sep = "\t")
if (nrow(carrier_cmp_dt) > 0) fwrite(carrier_cmp_dt, f_carrier_cmp, sep = "\t")

log_msg("Wrote diagnostics tables")

## Plot bookkeeping manifest
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
plot_manifest <- data.table(
    plot_file = character(),
    plot_id = character(),
    profile = character(),
    cohort_label = character(),
    focus_donor = character(),
    source_data = character(),
    description = character()
)
add_plot_manifest <- function(plot_file, plot_id, profile, source_data, description) {
    plot_manifest <<- rbind(
        plot_manifest,
        data.table(
            plot_file = plot_file,
            plot_id = plot_id,
            profile = profile,
            cohort_label = cohort_label,
            focus_donor = focus_donor,
            source_data = source_data,
            description = description
        ),
        fill = TRUE
    )
}

cohort_stamp <- sprintf("[%s / focus %s]", cohort_label, focus_donor)
make_full_plots <- plot_profile == "full"
focus_row <- diag_dt[BrNum == focus_donor]
if (nrow(focus_row) != 1L) stop("Focus donor not found in diagnostics table: ", focus_donor)

## Common plotting data for PCA-based panels
pc_dt <- data.table(
    BrNum = rownames(pc$x),
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    carrier = diag_dt$carrier[match(rownames(pc$x), diag_dt$BrNum)]
)
pc_dt <- merge(pc_dt, diag_dt[, .(BrNum, label_group)], by = "BrNum", all.x = TRUE, sort = FALSE)
pc_dt[, is_focus := BrNum == focus_donor]
pc_dt[, is_labeled := BrNum %in% label_donors]
pc_dt[, label_text := fifelse(BrNum == focus_donor, paste0(BrNum, " <- focus"), fifelse(is_labeled == TRUE, BrNum, NA_character_))]

if (make_full_plots) {
    ## 1) PCA plot
    p_pca <- ggplot(pc_dt, aes(x = PC1, y = PC2, color = carrier)) +
        geom_point(size = 2.5, alpha = 0.85) +
        geom_point(data = pc_dt[is_focus == TRUE], shape = 21, stroke = 1.1, color = "black", fill = "gold", size = 4) +
        geom_text(data = pc_dt[is_labeled == TRUE & is_focus == FALSE], aes(label = label_text), nudge_y = 0.45, size = 3.1, check_overlap = TRUE) +
        geom_text(data = pc_dt[is_focus == TRUE], aes(label = label_text), nudge_y = 0.6, fontface = "bold", size = 3.5) +
        theme_bw(base_size = 12) +
        labs(
            title = sprintf("%s PCA (21 donors)", cohort_stamp),
            subtitle = "logCPM of filterByExpr-passed genes",
            x = "PC1",
            y = "PC2",
            color = "Carrier"
        )
    f_pca <- file.path(plot_dir, sprintf("01_pca_focus_%s.png", focus_donor))
    ggsave(f_pca, p_pca, width = 7.4, height = 5.6, dpi = 140)
    add_plot_manifest(f_pca, "01", "full", "res21/logcpm21", "Global expression PCA with focus donor highlight.")

    ## 2) Median correlation rank plot
    corr_plot_dt <- copy(diag_dt[order(median_spearman_corr)])
    corr_plot_dt[, BrNum := factor(BrNum, levels = BrNum)]
    corr_plot_dt[, label_me := as.character(BrNum) %in% label_donors]
    p_corr <- ggplot(corr_plot_dt, aes(x = BrNum, y = median_spearman_corr, fill = label_group)) +
        geom_col() +
        geom_text(data = corr_plot_dt[label_me == TRUE], aes(label = as.character(BrNum)), hjust = -0.06, size = 2.9) +
        coord_flip() +
        scale_fill_manual(
            values = c("focus" = "#c0392b", "requested" = "#d35400", "outlier" = "#7d3c98", "near-outlier" = "#5dade2", "other" = "#7f8c8d"),
            guide = "none",
            drop = FALSE
        ) +
        expand_limits(y = max(corr_plot_dt$median_spearman_corr, na.rm = TRUE) + 0.02) +
        theme_bw(base_size = 11) +
        labs(
            title = sprintf("%s Median sample-sample Spearman correlation", cohort_stamp),
            x = "Donor",
            y = "Median correlation"
        )
    f_corr <- file.path(plot_dir, sprintf("02_median_corr_focus_%s.png", focus_donor))
    ggsave(f_corr, p_corr, width = 8.0, height = 6.8, dpi = 140)
    add_plot_manifest(f_corr, "02", "full", "res21/logcpm21", "Neighborhood similarity ranking across 21 donors.")

    ## 3) Focus z-distribution
    z_plot_dt <- rbindlist(list(
        data.table(group = "All tested genes", z = z_all),
        data.table(group = "20-donor DEGs (FDR<0.05)", z = z_deg20)
    ), fill = TRUE)
    p_z <- ggplot(z_plot_dt[is.finite(z)], aes(x = z, color = group, fill = group)) +
        geom_density(alpha = 0.2, linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        theme_bw(base_size = 12) +
        labs(
            title = sprintf("%s Focus donor directionality within %s", cohort_stamp, focus_carrier),
            subtitle = "Positive z means donor is higher than its own carrier-group mean",
            x = "Gene-wise z score",
            y = "Density",
            color = "Gene set",
            fill = "Gene set"
        )
    f_z <- file.path(plot_dir, sprintf("03_focus_z_density_%s.png", focus_donor))
    ggsave(f_z, p_z, width = 8.2, height = 5.4, dpi = 140)
    add_plot_manifest(f_z, "03", "full", "gene_dt/z summaries", "Focus-donor z density on all genes and DEG20 genes.")

    ## 4) Contrast gain vs DEG significance
    p_gain <- ggplot(
        gene_dt[is.finite(contrast_gain_drop_focus) & is.finite(FDR_20)],
        aes(x = contrast_gain_drop_focus, y = -log10(FDR_20 + 1e-300), color = is_deg20)
    ) +
        geom_point(alpha = 0.4, size = 0.9) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        theme_bw(base_size = 12) +
        labs(
            title = sprintf("%s Contrast gain when dropping %s", cohort_stamp, focus_donor),
            subtitle = "Each point is a gene; y-axis is 20-donor significance",
            x = "Carrier contrast gain (without focus donor minus with focus donor)",
            y = "-log10(FDR_20)",
            color = "DEG20"
        )
    f_gain <- file.path(plot_dir, sprintf("04_contrast_gain_vs_fdr20_%s.png", focus_donor))
    ggsave(f_gain, p_gain, width = 8.2, height = 5.4, dpi = 140)
    add_plot_manifest(f_gain, "04", "full", "gene_dt + res20", "Gene-level contrast gain versus significance in 20-donor model.")

    ## 5) Covariate z-score heatmap across donors
    covz_long <- melt(
        diag_dt[, .(BrNum, Age_z, Anc_Afr_z, pseudo_expr_chrM_ratio_z)],
        id.vars = "BrNum",
        variable.name = "covariate",
        value.name = "z"
    )
    covz_long[, covariate := factor(covariate, levels = c("Age_z", "Anc_Afr_z", "pseudo_expr_chrM_ratio_z"))]
    donor_order <- diag_dt[order(low_corr_rank, -pca_dist_pc123), BrNum]
    covz_long[, BrNum := factor(BrNum, levels = donor_order)]
    p_covz <- ggplot(covz_long, aes(x = covariate, y = BrNum, fill = z)) +
        geom_tile(color = "white", linewidth = 0.2) +
        scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0) +
        scale_y_discrete(labels = donor_axis_labels[donor_order]) +
        theme_bw(base_size = 11) +
        labs(
            title = sprintf("%s Covariate z-scores by donor", cohort_stamp),
            subtitle = "Donors ordered by expression outlier tendency",
            x = "Covariate",
            y = "Donor",
            fill = "z score"
        )
    f_covz <- file.path(plot_dir, sprintf("05_covariate_z_heatmap_%s.png", focus_donor))
    ggsave(f_covz, p_covz, width = 8.4, height = 6.8, dpi = 140)
    add_plot_manifest(f_covz, "05", "full", "diag_dt", "Covariate outlier pattern across 21 donors.")

    ## 6) Focus-donor outlier percentile across metrics
    metric_map <- data.table(
        metric = c("mahalanobis_rank_high", "leverage_rank_high", "pca_dist_rank_high",
                   "low_corr_rank", "extreme_gene_burden_rank_high", "chrM_ratio_rank_high"),
        label = c("Cov Mahalanobis", "Design leverage", "PCA distance",
                  "Low median corr", "Extreme-gene burden", "chrM ratio")
    )
    rank_cols <- metric_map$metric
    focus_metric <- data.table::melt(
        focus_row[, ..rank_cols],
        measure.vars = rank_cols,
        variable.name = "metric",
        value.name = "rank"
    )
    focus_metric <- merge(focus_metric, metric_map, by = "metric", all.x = TRUE)
    focus_metric[, outlier_percentile := 1 - ((rank - 1) / max(1, (nrow(diag_dt) - 1)))]
    focus_metric[, label := factor(label, levels = metric_map$label)]
    p_focus <- ggplot(focus_metric, aes(x = label, y = outlier_percentile)) +
        geom_col(fill = "#1f618d", width = 0.68) +
        geom_text(aes(label = sprintf("rank %d/%d", as.integer(rank), nrow(diag_dt))), hjust = -0.1, size = 3.0) +
        coord_flip() +
        ylim(0, 1.05) +
        theme_bw(base_size = 11) +
        labs(
            title = sprintf("%s Focus donor outlier percentile by metric", cohort_stamp),
            subtitle = "1.0 = most outlying donor on that metric",
            x = "Metric",
            y = "Outlier percentile"
        )
    f_focus_profile <- file.path(plot_dir, sprintf("06_focus_outlier_profile_%s.png", focus_donor))
    ggsave(f_focus_profile, p_focus, width = 8.6, height = 5.8, dpi = 140)
    add_plot_manifest(f_focus_profile, "06", "full", "diag_dt", "Focus donor rank profile across outlier metrics.")

    ## 7) Within-carrier comparison on 20-donor DEGs
    f_carrier_plot <- NA_character_
    if (nrow(carrier_cmp_dt) > 0) {
        cplot <- copy(carrier_cmp_dt)
        cplot[, BrNum := factor(BrNum, levels = BrNum[order(median_z_deg20)])]
        cplot[, is_focus := BrNum == focus_donor]
        cplot[, label_me := as.character(BrNum) %in% label_donors]
        p_carrier <- ggplot(cplot, aes(x = BrNum, y = median_z_deg20, fill = is_focus)) +
            geom_col() +
            geom_text(data = cplot[label_me == TRUE], aes(label = as.character(BrNum)), hjust = -0.07, size = 2.9) +
            coord_flip() +
            scale_fill_manual(values = c("FALSE" = "#7f8c8d", "TRUE" = "#b03a2e"), guide = "none") +
            expand_limits(y = max(cplot$median_z_deg20, na.rm = TRUE) + 0.2) +
            theme_bw(base_size = 11) +
            labs(
                title = sprintf("%s Within-%s donor comparison on 20-donor DEGs", cohort_stamp, focus_carrier),
                subtitle = "Gene-wise z scores computed within focus-carrier donors",
                x = "Donor",
                y = "Median z on DEG20 genes"
            )
        f_carrier_plot <- file.path(plot_dir, sprintf("07_carrier_deg20_medianz_%s.png", focus_donor))
        ggsave(f_carrier_plot, p_carrier, width = 8.4, height = 5.8, dpi = 140)
        add_plot_manifest(f_carrier_plot, "07", "full", "carrier_cmp_dt", "Focus donor compared to same-carrier peers on DEG20 genes.")
    }
}

## 8) Donor influence fingerprint (21 -> 20), based only on overlap leave-one runs
f_l1_scan <- here("code", "22_Mediation", "out-dge-test", "outlier_scan_leave1_21to20.tsv")
tt21_small <- res21$tt[, .(gene, logFC_21 = logFC)]
l1_influence_dt <- data.table()
if (file.exists(f_l1_scan)) {
    l1_scan <- fread(f_l1_scan)
    l1_scan <- l1_scan[drop1 %in% overlap_donors]
    if (nrow(l1_scan) > 0) {
        l1_scan[, donor_set_ok := vapply(seq_len(.N), function(i) {
            donor_vec <- sort(strsplit(donor_str[i], ",", fixed = TRUE)[[1]])
            setequal(donor_vec, setdiff(overlap_donors, drop1[i]))
        }, logical(1))]
        if (any(!l1_scan$donor_set_ok)) stop("Found leave-one rows not matching 21-overlap donor sets.")
        l1_influence_dt <- rbindlist(lapply(seq_len(nrow(l1_scan)), function(i) {
            rr <- l1_scan[i]
            if (!file.exists(rr$out_file)) return(NULL)
            tt20_i <- fread(rr$out_file)[, .(gene = outcome_gene_id, logFC_20 = logFC)]
            mm <- merge(tt21_small, tt20_i, by = "gene")
            gain_i <- mm$logFC_20 - mm$logFC_21
            data.table(
                drop1 = rr$drop1,
                drop1_apoe = rr$drop1_apoe,
                drop1_carrier = rr$drop1_carrier,
                n_deg_fdr05 = rr$n_deg_fdr05,
                min_fdr = rr$min_fdr,
                n_genes_compared = nrow(mm),
                median_gain = median(gain_i, na.rm = TRUE),
                frac_gain_positive = mean(gain_i > 0, na.rm = TRUE),
                mean_abs_gain = mean(abs(gain_i), na.rm = TRUE)
            )
        }), fill = TRUE)
    }
}
if (nrow(l1_influence_dt) < 1) {
    log_msg("Leave-one summary file unavailable/incomplete; computing influence summaries directly.")
    l1_influence_dt <- rbindlist(lapply(overlap_donors, function(d) {
        rr <- run_model(setdiff(overlap_donors, d))
        mm <- merge(tt21_small, rr$tt[, .(gene, logFC_20 = logFC)], by = "gene")
        gain_i <- mm$logFC_20 - mm$logFC_21
        dm <- diag_dt[BrNum == d]
        data.table(
            drop1 = d,
            drop1_apoe = dm$APOE_syn[1],
            drop1_carrier = dm$carrier[1],
            n_deg_fdr05 = sum(rr$tt$adj.P.Val < 0.05, na.rm = TRUE),
            min_fdr = min(rr$tt$adj.P.Val, na.rm = TRUE),
            n_genes_compared = nrow(mm),
            median_gain = median(gain_i, na.rm = TRUE),
            frac_gain_positive = mean(gain_i > 0, na.rm = TRUE),
            mean_abs_gain = mean(abs(gain_i), na.rm = TRUE)
        )
    }), fill = TRUE)
}
setorder(l1_influence_dt, -median_gain, -frac_gain_positive)
fwrite(l1_influence_dt, f_l1_influence, sep = "\t")

l1_influence_dt[, is_focus := drop1 == focus_donor]
l1_influence_dt[, label_me := (drop1 %in% label_donors) | (n_deg_fdr05 > 0)]
p_influence <- ggplot(l1_influence_dt, aes(x = median_gain, y = frac_gain_positive, color = n_deg_fdr05, shape = drop1_carrier)) +
    geom_point(size = 3.0, alpha = 0.9) +
    geom_point(data = l1_influence_dt[is_focus == TRUE], shape = 21, fill = "gold", color = "black", size = 4.3, stroke = 1.1) +
    geom_text(data = l1_influence_dt[label_me == TRUE], aes(label = drop1), nudge_y = 0.01, size = 3.0, check_overlap = TRUE) +
    scale_color_gradient(low = "#5dade2", high = "#b03a2e") +
    theme_bw(base_size = 12) +
    labs(
        title = sprintf("%s Donor influence fingerprint (21->20)", cohort_stamp),
        subtitle = "Each point summarizes gene-level carrier-contrast gain for one dropped donor",
        x = "Median carrier-contrast gain across genes",
        y = "Fraction genes with positive gain",
        color = "DEG count",
        shape = "Dropped donor carrier"
    )
f_influence_plot <- file.path(plot_dir, "08_donor_influence_fingerprint_21to20.png")
ggsave(f_influence_plot, p_influence, width = 8.4, height = 5.8, dpi = 140)
add_plot_manifest(f_influence_plot, "08", "explanatory", "outlier_scan_leave1_21to20 + res21", "Donor-wise directional attenuation summary in leave-one-out runs.")

## 9) DEG20 directionality heatmap for all 21 donors
rd_map <- data.table(
    gene = rownames(sce_base),
    gene_name = as.character(SummarizedExperiment::rowData(sce_base)$gene_name),
    gene_type = as.character(SummarizedExperiment::rowData(sce_base)$gene_type)
)
top_heat_genes <- gene_dt[is_deg20 == TRUE & is.finite(contrast_gain_drop_focus)][order(-contrast_gain_drop_focus)][1:min(40, .N), gene]
if (length(top_heat_genes) < 1) {
    top_heat_genes <- gene_dt[is.finite(contrast_gain_drop_focus)][order(-contrast_gain_drop_focus)][1:min(40, .N), gene]
}
carrier_by_col <- diag_dt$carrier[match(colnames(logcpm21), diag_dt$BrNum)]
zmat_carrier_all <- matrix(NA_real_, nrow = nrow(logcpm21), ncol = ncol(logcpm21), dimnames = dimnames(logcpm21))
for (cg in unique(carrier_by_col)) {
    cidx <- which(carrier_by_col == cg)
    subm <- logcpm21[, cidx, drop = FALSE]
    mu <- rowMeans(subm, na.rm = TRUE)
    sig <- apply(subm, 1, sd, na.rm = TRUE)
    sig[!is.finite(sig) | sig == 0] <- NA_real_
    zsub <- sweep(subm, 1, mu, "-")
    zsub <- sweep(zsub, 1, sig, "/")
    zmat_carrier_all[, cidx] <- zsub
}
donor_levels <- diag_dt[order(carrier, low_corr_rank), BrNum]
donor_labels <- donor_axis_labels[donor_levels]
gene_label_dt <- merge(data.table(gene = top_heat_genes), rd_map[, .(gene, gene_name)], by = "gene", all.x = TRUE)
gene_label_dt[, gene_label := ifelse(is.na(gene_name) | gene_name == "", gene, paste0(gene_name, " (", gene, ")"))]
heat_dt <- as.data.table(as.table(zmat_carrier_all[top_heat_genes, donor_levels, drop = FALSE]))
setnames(heat_dt, c("gene", "BrNum", "z"))
heat_dt <- merge(heat_dt, diag_dt[, .(BrNum, carrier)], by = "BrNum", all.x = TRUE)
heat_dt <- merge(heat_dt, gene_label_dt[, .(gene, gene_label)], by = "gene", all.x = TRUE)
heat_dt[, gene_label := factor(gene_label, levels = rev(gene_label_dt$gene_label))]
heat_dt[, BrNum := factor(BrNum, levels = donor_levels)]
p_heat <- ggplot(heat_dt, aes(x = gene_label, y = BrNum, fill = z)) +
    geom_tile() +
    facet_grid(carrier ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0, na.value = "grey92") +
    scale_y_discrete(labels = setNames(donor_labels, donor_levels)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
    labs(
        title = sprintf("%s DEG20 directionality heatmap (carrier-normalized z)", cohort_stamp),
        subtitle = "Top genes ranked by contrast gain after dropping focus donor",
        x = "Genes",
        y = "Donor",
        fill = "z score"
    )
f_heat_plot <- file.path(plot_dir, sprintf("09_deg20_directionality_heatmap_%s.png", focus_donor))
ggsave(f_heat_plot, p_heat, width = 11.2, height = 6.8, dpi = 140)
add_plot_manifest(f_heat_plot, "09", "explanatory", "logcpm21 + DEG20 genes", "Gene-by-donor carrier-normalized directionality on top DEG20 genes.")

## 10) Carrier-aware expression positioning map
p_pca_ellipse <- ggplot(pc_dt, aes(x = PC1, y = PC2, color = carrier)) +
    stat_ellipse(aes(group = carrier), type = "norm", level = 0.80, linewidth = 0.8, linetype = 2, show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_point(data = pc_dt[is_focus == TRUE], shape = 21, stroke = 1.2, color = "black", fill = "gold", size = 4.4) +
    geom_text(data = pc_dt[is_labeled == TRUE & is_focus == FALSE], aes(label = label_text), nudge_y = 0.45, size = 3.0, check_overlap = TRUE) +
    geom_text(data = pc_dt[is_focus == TRUE], aes(label = label_text), nudge_y = 0.7, fontface = "bold", size = 3.5) +
    theme_bw(base_size = 12) +
    labs(
        title = sprintf("%s Carrier-aware expression positioning", cohort_stamp),
        subtitle = "PCA with carrier ellipses (80%) and focus donor highlight",
        x = "PC1",
        y = "PC2",
        color = "Carrier"
    )
f_pca_ellipse <- file.path(plot_dir, sprintf("10_pca_ellipse_carrier_%s.png", focus_donor))
ggsave(f_pca_ellipse, p_pca_ellipse, width = 8.0, height = 5.8, dpi = 140)
add_plot_manifest(f_pca_ellipse, "10", "explanatory", "res21/logcpm21", "PCA map with carrier ellipses and focus donor highlight.")

## 11) Biology coherence via hallmark enrichment on gain-ranked genes
fgsea_dt <- data.table()
rank_dt <- merge(gene_dt[, .(gene, contrast_gain_drop_focus)], rd_map[, .(gene, gene_name)], by = "gene", all.x = TRUE)
rank_dt <- rank_dt[is.finite(contrast_gain_drop_focus) & !is.na(gene_name) & gene_name != ""]
rank_sym <- rank_dt[, .(stat = mean(contrast_gain_drop_focus, na.rm = TRUE)), by = gene_name]
rank_sym <- rank_sym[is.finite(stat)]
if (nrow(rank_sym) > 50) {
    rank_vec <- rank_sym$stat
    names(rank_vec) <- rank_sym$gene_name
    rank_vec <- sort(rank_vec, decreasing = TRUE)
    msig_h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    pathways_h <- split(msig_h$gene_symbol, msig_h$gs_name)
    fg <- fgsea::fgsea(pathways = pathways_h, stats = rank_vec, minSize = 10, maxSize = 500, nPermSimple = 10000)
    fgsea_dt <- as.data.table(fg)[order(padj, -abs(NES))]
}
if (nrow(fgsea_dt) > 0) fwrite(fgsea_dt, f_fgsea, sep = "\t")

if (nrow(fgsea_dt) > 0) {
    fg_plot <- fgsea_dt[is.finite(NES) & is.finite(padj)][order(padj, -abs(NES))][1:min(15, .N)]
    fg_plot[, pathway := factor(pathway, levels = rev(pathway))]
    p_enrich <- ggplot(fg_plot, aes(x = NES, y = pathway, size = -log10(padj + 1e-300), color = padj)) +
        geom_point(alpha = 0.9) +
        scale_color_gradient(low = "#b2182b", high = "#2166ac", trans = "reverse") +
        theme_bw(base_size = 11) +
        labs(
            title = sprintf("%s High-gain biology coherence (Hallmark FGSEA)", cohort_stamp),
            subtitle = "Pathways enriched among genes ranked by contrast gain after dropping focus donor",
            x = "Normalized enrichment score (NES)",
            y = "Hallmark pathway",
            size = "-log10(FDR)",
            color = "FDR"
        )
} else {
    p_enrich <- ggplot(data.table(x = 1, y = 1), aes(x, y)) +
        geom_text(label = "FGSEA enrichment unavailable or not significant", size = 4.2) +
        theme_void() +
        labs(title = sprintf("%s High-gain biology coherence", cohort_stamp))
}
f_enrich_plot <- file.path(plot_dir, sprintf("11_high_gain_enrichment_dotplot_%s.png", focus_donor))
ggsave(f_enrich_plot, p_enrich, width = 9.2, height = 5.8, dpi = 140)
add_plot_manifest(f_enrich_plot, "11", "explanatory", "gene_dt + msigdbr(H) + fgsea", "Hallmark enrichment on genes ranked by contrast gain.")

fwrite(plot_manifest, f_plot_manifest, sep = "\t")
log_msg("Saved plots and plot manifest")

## Scope readme for output interpretation
scope_txt <- c(
    sprintf("Scope map generated: %s", Sys.time()),
    "",
    "1) Stepdown diagnostics (expanded cohort -> overlap cohort):",
    "  - Source: baseline_testing.R run_mode=stepdown",
    "  - Typical location: out-dge-test/plots/01..07 and run10_21donor_diagnostic_report.txt",
    "",
    "2) Outlier/sweep scans (21->20/19):",
    "  - Source: baseline_testing.R run_mode=outlier or sweep",
    "  - Typical location: out-dge-test/outlier_scan_*.tsv and plots/10,11",
    "",
    sprintf("3) %s focused EDA (21->20 only, focus donor):", focus_donor),
    "  - Source: eda_check.R",
    sprintf("  - Location: %s", out_dir),
    "  - This folder should contain only 21-overlap/20-drop-focus artifacts and no stepdown-cohort claims."
)
writeLines(scope_txt, con = f_scope)

## Text report with explicit cohort metadata
report_file <- file.path(out_dir, sprintf("eda_report_%s.txt", focus_donor))
con <- file(report_file, "wt")
writeLines(sprintf("EDA report for donor impact check (%s)", Sys.time()), con)
writeLines(sprintf("Focus donor: %s", focus_donor), con)
writeLines(sprintf("Cohort label: %s", cohort_label), con)
writeLines(sprintf("Overlap donors observed: %d ; asserted expected: %d", length(overlap_donors), assert_n_overlap), con)
writeLines(sprintf("Overlap donor hash: %s", overlap_donor_hash), con)
writeLines(sprintf("Plot profile: %s", plot_profile), con)
writeLines("Rank convention: for outlier ranks, 1/N = most outlying and N/N = least outlying.", con)
writeLines("", con)
writeLines("1) Model rerun summary", con)
writeLines(sprintf("  21 overlap donors: %d ; 20 donors after dropping %s: %d", length(overlap_donors), focus_donor, length(donors20)), con)
writeLines(sprintf("  21-donor DEGs (FDR<0.05): %d ; min FDR: %s", n_deg21, fmt(min_fdr21, 6)), con)
writeLines(sprintf("  20-donor DEGs (FDR<0.05): %d ; min FDR: %s", n_deg20, fmt(min_fdr20, 6)), con)
writeLines("", con)

writeLines("2) Focus donor diagnostics in 21-donor set", con)
writeLines(sprintf("  APOE=%s ; carrier=%s ; Sex=%s ; exp_round=%s",
                   focus_row$APOE_syn[1], focus_row$carrier[1], focus_row$Sex[1], focus_row$exp_round[1]), con)
writeLines(sprintf("  Age=%.2f (z=%s), Anc_Afr=%.4f (z=%s), chrM_ratio=%.6f (z=%s; rank-high=%d/%d)",
                   focus_row$Age[1], fmt(focus_row$Age_z[1]), focus_row$Anc_Afr[1], fmt(focus_row$Anc_Afr_z[1]),
                   focus_row$pseudo_expr_chrM_ratio[1], fmt(focus_row$pseudo_expr_chrM_ratio_z[1]),
                   focus_row$chrM_ratio_rank_high[1], nrow(diag_dt)), con)
writeLines(sprintf("  Covariate Mahalanobis rank (high=outlying): %d/%d",
                   focus_row$mahalanobis_rank_high[1], nrow(diag_dt)), con)
writeLines(sprintf("  Design leverage rank (high=outlying): %d/%d",
                   focus_row$leverage_rank_high[1], nrow(diag_dt)), con)
writeLines(sprintf("  PCA distance rank (high=outlying): %d/%d",
                   focus_row$pca_dist_rank_high[1], nrow(diag_dt)), con)
writeLines(sprintf("  Median Spearman correlation rank (low=outlying): %d/%d",
                   focus_row$low_corr_rank[1], nrow(diag_dt)), con)
writeLines(sprintf("  Extreme-gene burden rank (high=outlying): %d/%d",
                   focus_row$extreme_gene_burden_rank_high[1], nrow(diag_dt)), con)
writeLines("", con)

writeLines("3) Gene-level pattern linked to attenuation", con)
writeLines(sprintf("  Focus donor carrier-group (%s) z-score across all genes:", focus_carrier), con)
writeLines(sprintf("    median=%s ; mean=%s ; frac(z>0)=%s ; frac(z>1.5)=%s ; frac(z<-1.5)=%s",
                   fmt(median(z_all, na.rm = TRUE)), fmt(mean(z_all, na.rm = TRUE)),
                   fmt(mean(z_all > 0, na.rm = TRUE), 4), fmt(mean(z_all > 1.5, na.rm = TRUE), 4),
                   fmt(mean(z_all < -1.5, na.rm = TRUE), 4)), con)
writeLines(sprintf("    binomial test p(z>0 > 0.5): p=%s", format(sign_test_all$p.value, scientific = TRUE, digits = 4)), con)
if (length(z_deg20) > 0) {
    writeLines(sprintf("  In 20-donor DEGs only (n=%d):", length(z_deg20)), con)
    writeLines(sprintf("    median z=%s ; mean z=%s ; frac(z>0)=%s ; frac(z>1.5)=%s",
                       fmt(median(z_deg20, na.rm = TRUE)), fmt(mean(z_deg20, na.rm = TRUE)),
                       fmt(mean(z_deg20 > 0, na.rm = TRUE), 4), fmt(mean(z_deg20 > 1.5, na.rm = TRUE), 4)), con)
    writeLines(sprintf("    median contrast gain when dropping focus donor=%s ; frac(gain>0)=%s",
                       fmt(median(gain_deg20, na.rm = TRUE), 4), fmt(mean(gain_deg20 > 0, na.rm = TRUE), 4)), con)
}
writeLines("", con)

writeLines("4) Interpretation notes", con)
writeLines("  - Observation (21->20 DEG jump) is separated from explanatory analyses.", con)
writeLines("  - Expression-position and directional attenuation are the primary explanatory layers.", con)
writeLines("  - Mahalanobis and covariates are treated as rule-out checks, not primary evidence.", con)
writeLines("  - This report is 21-overlap focused; stepdown-cohort artifacts are outside this folder.", con)
writeLines("", con)

writeLines("5) Output files", con)
writeLines(sprintf("  %s", f_diag), con)
writeLines(sprintf("  %s", f_focus), con)
writeLines(sprintf("  %s", f_gene), con)
writeLines(sprintf("  %s", f_sum), con)
writeLines(sprintf("  %s", f_top_deg20), con)
writeLines(sprintf("  %s", f_top_gain), con)
if (nrow(carrier_cmp_dt) > 0) writeLines(sprintf("  %s", f_carrier_cmp), con)
writeLines(sprintf("  %s", f_l1_influence), con)
if (nrow(fgsea_dt) > 0) writeLines(sprintf("  %s", f_fgsea), con)
writeLines(sprintf("  %s", f_plot_manifest), con)
writeLines(sprintf("  %s", f_scope), con)
writeLines(sprintf("  %s", log_file), con)
close(con)

## Digest report (ordered narrative with rule-out placement)
e2_rank_txt <- "NA"
if (nrow(carrier_cmp_dt) > 0 && focus_donor %in% carrier_cmp_dt$BrNum) {
    fr <- carrier_cmp_dt[BrNum == focus_donor]
    e2_rank_txt <- sprintf("%d/%d", fr$median_z_rank_high[1], nrow(carrier_cmp_dt))
}
rank_desc <- function(rank_val, n_total) {
    if (!is.finite(rank_val) || !is.finite(n_total) || n_total < 2) return("not estimable")
    if (rank_val <= 3) return("strong outlier")
    if (rank_val <= 5) return("near outlier")
    if (rank_val <= 10) return("mid-range")
    "not outlying"
}
frac_pos_gain <- if (length(gain_deg20) > 0) mean(gain_deg20 > 0, na.rm = TRUE) else NA_real_
median_gain_deg20 <- if (length(gain_deg20) > 0) median(gain_deg20, na.rm = TRUE) else NA_real_
median_z_deg20 <- if (length(z_deg20) > 0) median(z_deg20, na.rm = TRUE) else NA_real_
strong_expr_outlier <- (focus_row$pca_dist_rank_high[1] <= 5L) || (focus_row$low_corr_rank[1] <= 5L) || (focus_row$extreme_gene_burden_rank_high[1] <= 5L)
directional_positive <- is.finite(median_gain_deg20) && is.finite(frac_pos_gain) && (median_gain_deg20 > 0) && (frac_pos_gain > 0.5)
directional_negative <- is.finite(median_gain_deg20) && is.finite(frac_pos_gain) && (median_gain_deg20 < 0) && (frac_pos_gain < 0.5)
integrated_conclusion <- if (n_deg20 < 1) {
    sprintf("Dropping %s does not create DEGs, so there is no practical donor-driven attenuation signal to explain.", focus_donor)
} else if (directional_positive && strong_expr_outlier) {
    sprintf("%s shows a coherent expression-space outlier pattern and a donor-direction shift consistent with attenuation of the carrier contrast.", focus_donor)
} else if (directional_negative && !strong_expr_outlier) {
    sprintf("%s is not a strong expression-space outlier, and the 21->20 DEG emergence appears localized/threshold-based rather than broad directional attenuation.", focus_donor)
} else if (directional_negative && strong_expr_outlier) {
    sprintf("%s has some expression outlier characteristics, but the emergent genes show mixed or opposite-direction shifts; treat exclusion as sensitivity analysis, not a primary mechanistic rule.", focus_donor)
} else {
    sprintf("%s has a mixed donor-impact signature: DEG emergence occurs, but directional attenuation evidence is not uniformly strong.", focus_donor)
}
decision_statement <- if (n_deg20 < 1) {
    sprintf("No sensitivity exclusion case for %s: dropping this donor does not recover DEGs.", focus_donor)
} else if (directional_positive && strong_expr_outlier) {
    sprintf("Sensitivity exclusion of %s is mechanistically supported; not based on DEG count alone.", focus_donor)
} else {
    sprintf("For %s, any exclusion should be sensitivity-only: DEG emergence is present, but broad directional attenuation evidence is limited.", focus_donor)
}
mechanism_line <- if (n_deg20 < 1) {
    "no emergent DEGs after donor removal."
} else if (directional_positive) {
    sprintf("emergent genes show positive directional gain (median gain=%s; frac(gain>0)=%s).", fmt(median_gain_deg20, 3), fmt(frac_pos_gain, 3))
} else if (directional_negative) {
    sprintf("emergent genes do not show positive directional attenuation overall (median gain=%s; frac(gain>0)=%s).", fmt(median_gain_deg20, 3), fmt(frac_pos_gain, 3))
} else {
    sprintf("emergent genes show mixed directionality (median gain=%s; frac(gain>0)=%s).", fmt(median_gain_deg20, 3), fmt(frac_pos_gain, 3))
}
digest_lines <- c(
    sprintf("# %s donor-impact analysis: coherent plain-language narrative", focus_donor),
    "",
    "## Core question",
    sprintf("Why does removing %s change the APOE carrier analysis from no DEGs to many DEGs?", focus_donor),
    "",
    "## Step 0: Observation (not explanation)",
    sprintf("- 21 donors: **%d DEGs** (best FDR = **%s**)", n_deg21, fmt(min_fdr21, 4)),
    sprintf("- Drop %s (20 donors): **%d DEGs** (best FDR = **%s**)", focus_donor, n_deg20, fmt(min_fdr20, 5)),
    "",
    "## Step 1: Is the focus donor globally unusual in expression space?",
    sprintf("- PCA distance rank: **%d/%d**", focus_row$pca_dist_rank_high[1], nrow(diag_dt)),
    sprintf("- Low-correlation rank: **%d/%d**", focus_row$low_corr_rank[1], nrow(diag_dt)),
    sprintf("- Extreme-gene burden rank: **%d/%d**", focus_row$extreme_gene_burden_rank_high[1], nrow(diag_dt)),
    "",
    "## Step 2: Could this mainly be a simple covariate artifact?",
    "- Rank convention reminder: **1/N means most outlying** on that metric.",
    sprintf("- Mahalanobis rank: **%d/%d** (%s)", focus_row$mahalanobis_rank_high[1], nrow(diag_dt), rank_desc(focus_row$mahalanobis_rank_high[1], nrow(diag_dt))),
    sprintf("- Leverage rank: **%d/%d** (%s)", focus_row$leverage_rank_high[1], nrow(diag_dt), rank_desc(focus_row$leverage_rank_high[1], nrow(diag_dt))),
    sprintf("- chrM ratio rank: **%d/%d** (%s)", focus_row$chrM_ratio_rank_high[1], nrow(diag_dt), rank_desc(focus_row$chrM_ratio_rank_high[1], nrow(diag_dt))),
    "",
    "## Step 3: Directional shift on emergent DEGs",
    sprintf("- Fraction with positive contrast gain after dropping %s: **%s**", focus_donor, fmt(frac_pos_gain, 3)),
    sprintf("- Median contrast gain: **%s**", fmt(median_gain_deg20, 3)),
    sprintf("- Focus donor DEG20 median z (within carrier): **%s**", fmt(median_z_deg20, 3)),
    "",
    "## Step 4: Same-carrier context",
    sprintf("- Focus donor median-z rank on DEG20 genes within %s donors: **%s**", focus_carrier, e2_rank_txt),
    "",
    "## Integrated conclusion",
    integrated_conclusion
)
writeLines(digest_lines, con = f_digest)

## One-slide markdown artifact
slide_lines <- c(
    sprintf("# Why %s suppresses APOE-carrier DE signal (21-overlap analysis)", focus_donor),
    "",
    sprintf("- **Observation only:** 21 donors -> 20 donors (drop %s) changes DEGs from **%d** to **%d**.", focus_donor, n_deg21, n_deg20),
    sprintf("- **Mechanism evidence:** %s", mechanism_line),
    sprintf("- **Positioning evidence:** PCA distance rank=%d/%d and low-correlation rank=%d/%d.", focus_row$pca_dist_rank_high[1], nrow(diag_dt), focus_row$low_corr_rank[1], nrow(diag_dt)),
    sprintf("- **Rule-out check:** Mahalanobis rank=%d/%d (%s; rank 1 = most outlying).", focus_row$mahalanobis_rank_high[1], nrow(diag_dt), rank_desc(focus_row$mahalanobis_rank_high[1], nrow(diag_dt))),
    sprintf("- **Same-carrier context:** %s median-z rank on DEG20 genes within %s donors = %s.", focus_donor, focus_carrier, e2_rank_txt),
    "",
    sprintf("**Decision statement:** %s", decision_statement),
    "",
    "Recommended visuals:",
    sprintf("1. `%s`", f_influence_plot),
    sprintf("2. `%s`", f_pca_ellipse),
    sprintf("3. `%s`", f_heat_plot)
)
writeLines(slide_lines, con = f_slide_md)
slide_html <- file.path(out_dir, sprintf("%s_one_slide.html", tolower(focus_donor)))
slide_pptx <- file.path(out_dir, sprintf("%s_one_slide.pptx", tolower(focus_donor)))
if (nzchar(Sys.which("pandoc"))) {
    try(suppressWarnings(system2("pandoc", args = c(f_slide_md, "-s", "-o", slide_html))), silent = TRUE)
    try(suppressWarnings(system2("pandoc", args = c(f_slide_md, "-t", "pptx", "-o", slide_pptx))), silent = TRUE)
}

log_msg("Wrote report: ", report_file)
log_msg("Wrote digest: ", f_digest)
log_msg("Wrote slide markdown: ", f_slide_md)
log_msg("Done")
