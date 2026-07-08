#' Plot outcome gene expression with and without adjusting for a mediator gene
#'
#' Adapted from `plot_DEG_express()` for the Baron & Kenny mediation
#' workflow. Where `plot_DEG_express()` shows one cleaned expression panel
#' per outcome gene (covariates regressed out), this function shows TWO
#' panels side by side for each outcome gene (or THREE, if
#' `plot_mediator_panel = TRUE`):
#'   - (optional) "mediator" panel: the mediator gene itself, in its own
#'     cluster, plotted with `plot_DEG_express()` -- this is the Step 2
#'     (X->M) view (e.g. NPTXR in Astro.1).
#'   - "unadjusted": outcome gene expression with the usual nuisance
#'     covariates (Sex/Age/Anc_Afr/etc, per `mod`/`cleanY_P`) regressed out
#'     -- this is the Step 1 (X->Y) view, identical to `plot_DEG_express()`.
#'   - "adjusted for {mediator}": the same, but with the mediator gene's
#'     expression ALSO regressed out as an additional covariate -- this is
#'     the Step 3b (X->Y | M) view. If the carrier-vs-outcome relationship
#'     visibly flattens between the two panels, that's the attenuation
#'     signature of mediation.
#'
#' The mediator's expression comes from a second (typically different
#' cluster) pseudobulk object, matched to the outcome object by BrNum.
#' Only donors present in both objects are plotted, so both panels show
#' the exact same sample set.
#'
#' @param sce [SingleCellExperiment] outcome pseudobulk object.
#' @param sce_mediator [SingleCellExperiment] mediator pseudobulk object
#' (can be the same object as `sce` if mediator and outcome share a
#' pseudobulk file, e.g. the `Xenium_Oligo.3_Astro` scenarios).
#' @param stats A `data.frame()` with mediation validation results -- one
#' row per outcome gene for the given `mediator`/`med_clus`, matching the
#' shape of `mediation_validation_details()`'s output / the
#' `mediation_summary-all_scenarios_*.csv` files (columns `outcome`,
#' `mediator`, `P.Value_Xen`, `t_Xen`, `P.Value_Xen_carrier`,
#' `t_Xen_carrier`, `P.Value_Xen_med`, `t_Xen_med`, `Xen_mediated`).
#' @param clus A `character(1)` outcome cluster (matches `cluster_col` in
#' `colData(sce)`).
#' @param med_clus A `character(1)` mediator cluster (matches
#' `med_cluster_col` in `colData(sce_mediator)`). This is the cluster
#' label actually used to pull mediator expression, e.g. `med_cl_test`
#' in the APOE-stratified scenarios.
#' @param mediator_gene A `character(1)` mediator gene name (must be in
#' `rownames(sce_mediator)` after `gene_col` is applied, and must match
#' `stats$mediator`).
#' @param gene A `character()` vector of outcome gene(s) to plot (subset
#' of `stats$outcome`).
#' @param outcome_col,mediator_col The `character(1)` column names in
#' `stats` identifying the outcome gene and mediator gene respectively.
#' @param gene_col The `character(1)` name of the `rowData()` column
#' (same on both `sce` and `sce_mediator`) giving the gene symbol used to
#' match against `gene`/`mediator_gene`.
#' @param cluster_col,med_cluster_col The `character(1)` `colData()`
#' column names identifying cluster membership on `sce`/`sce_mediator`
#' respectively.
#' @param category_col The `character(1)` `colData()` column to color/
#' group boxplots by (e.g. `"APOE_carrier"`).
#' @param mod A `formula` or design `matrix` used to clean the OUTCOME
#' gene's expression, same convention as `plot_DEG_express()` -- the
#' first `cleanY_P` columns are kept (not regressed out); everything
#' after is treated as a nuisance covariate to remove. The mediator
#' column is appended as an additional trailing covariate for the
#' "adjusted" panel only, so `cleanY_P` does not need to change.
#' @param cleanY_P An `integer(1)`, number of leading model matrix
#' columns to keep (matches the covariate-of-interest columns, e.g. the
#' APOE genotype factor levels).
#' @param color_pal A named `character()` color palette matching
#' `category_col` values.
#' @param plot_points,ncol Passed through to `plot_gene_express()`.
#' @param min_donors An `integer(1)`, minimum number of donors shared
#' between `sce` and `sce_mediator` required to proceed (default 10,
#' matching the mediation fitting loop's own donor-count floor).
#' @param plot_mediator_panel A `logical(1)`. If `TRUE`, prepend a first
#' panel showing the mediator gene itself (via `plot_DEG_express()`) in
#' its own cluster, e.g. NPTXR in Astro.1 -- the Step 2 (X->M) view.
#' Requires `mediator_stats`.
#' @param mediator_stats A `data.frame()` in the standard
#' `plot_DEG_express()` stats shape (columns `cluster`, plus whatever
#' `mediator_pval_col`/`mediator_fc_col`/`gene_col` point to) -- the
#' mediator's own carrier-DE result within `med_clus`. Only used when
#' `plot_mediator_panel = TRUE`.
#' @param mediator_pval_col,mediator_fc_col Column names in
#' `mediator_stats` for the mediator panel's p-value/logFC annotation,
#' passed straight through to `plot_DEG_express()`.
#' @param mediator_mod,mediator_cleanY_P Optional `formula`/`matrix` and
#' `integer(1)` used to clean the mediator's own expression for its
#' panel. Default to `mod`/`cleanY_P` (i.e. assume the same design) if
#' left `NULL` -- override if the mediator's cluster needs a different
#' model.
#' @param anno_stat_suffix A `character(1)` identifying which stat
#' columns to read from `stats` for annotation, e.g. `"Xen"` (default)
#' for Xenium validation columns (`P.Value_Xen`, `P.Value_Xen_carrier`,
#' `P.Value_Xen_med`, `Xen_mediated`) or `"SN"` for discovery-side
#' columns (`fdr_SN`, `fdr_SN_carrier`, `fdr_SN_med`, `SN_mediated`).
#' @param signif_stat A `character(1)`, the p-value-like column prefix to
#' look for (paired with `anno_stat_suffix`), e.g. `"P.Value"` or
#' `"fdr"`.
#' @param signif_label A `character(1)` display label for the annotation
#' text (e.g. `"P"` vs `"fdr"`). Defaults to `"P"` when
#' `signif_stat == "P.Value"`, otherwise to `signif_stat` itself.
#'
#' Each of the three annotation pieces (Step 1 base, Step 3b
#' carrier-adjusted, Step 3a mediator-association) is only added if its
#' corresponding columns (`{signif_stat}_{anno_stat_suffix}[_carrier|_med]`,
#' `t_{anno_stat_suffix}[_carrier|_med]`) are present in `stats` -- e.g.
#' if Step 3b hasn't been computed yet, omit the `*_carrier` columns and
#' that piece of the annotation is silently skipped (with a `message()`)
#' rather than erroring.
#'
#' @return A `patchwork` object: two (or three, with the mediator panel)
#' `ggplot2` panels side by side, one facet per outcome gene within each
#' outcome panel.
#' @export
#'
#' @examples
#' ## Plot NPTXR (outcome, Oligo.3) with and without adjusting for
#' ## Astro.1's NPTXR-associated mediator gene
#' plot_DEG_mediated_express(
#'     sce = sce_pb,
#'     sce_mediator = sce_mediator_pb,
#'     stats = mediation_validation_details_test,
#'     clus = "Oligo.3",
#'     med_clus = "Astro.1",
#'     mediator_gene = "NPTXR",
#'     gene = c("ENC1", "GAD1", "SLC17A7"),
#'     gene_col = "gene_name",
#'     cluster_col = "registration_variable",
#'     med_cluster_col = "registration_variable",
#'     category_col = "APOE_carrier",
#'     mod = ~APOE_carrier_syn + Age + Anc_Afr,
#'     cleanY_P = 2
#' )
#' @family expression plotting functions
#' @importFrom ggplot2 ggplot geom_label labs
plot_DEG_mediated_express <- function(sce,
                                       sce_mediator,
                                       stats,
                                       clus,
                                       med_clus,
                                       mediator_gene,
                                       gene,
                                       outcome_col = "outcome",
                                       mediator_col = "mediator",
                                       gene_col = "gene_name",
                                       cluster_col = "registration_variable",
                                       med_cluster_col = "registration_variable",
                                       category_col = "APOE_carrier",
                                       mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
                                       cleanY_P = 4,
                                       color_pal = NULL,
                                       plot_points = FALSE,
                                       ncol = 2,
                                       min_donors = 10,
                                       plot_mediator_panel = FALSE,
                                       mediator_stats = NULL,
                                       mediator_pval_col = "vlmf_adj.P.Val",
                                       mediator_fc_col = "vlmf_logFC",
                                       mediator_mod = NULL,
                                       mediator_cleanY_P = NULL,
                                       anno_stat_suffix = "Xen",
                                       signif_stat = "P.Value",
                                       signif_label = NULL) {

    stopifnot(cluster_col %in% colnames(colData(sce)))
    stopifnot(med_cluster_col %in% colnames(colData(sce_mediator)))
    stopifnot(clus %in% sce[[cluster_col]])
    stopifnot(med_clus %in% sce_mediator[[med_cluster_col]])
    stopifnot(outcome_col %in% colnames(stats))
    stopifnot(mediator_col %in% colnames(stats))
    stopifnot(all(gene %in% stats[[outcome_col]]))
    stopifnot(mediator_gene %in% stats[[mediator_col]])

    # RCMD fix
    rank_int <- Symbol <- anno_str <- Var1 <- NULL

    ## ------------------------------------------------------------------
    ## Column-name pattern for annotation stats. Built from
    ## anno_stat_suffix/signif_stat so the same function works whether
    ## `stats` uses the Xenium validation naming (P.Value_Xen,
    ## P.Value_Xen_carrier, P.Value_Xen_med, Xen_mediated) or a
    ## discovery-side naming (fdr_SN, fdr_SN_carrier, fdr_SN_med,
    ## SN_mediated) -- pass anno_stat_suffix="SN", signif_stat="fdr" for
    ## the latter.
    ##
    ## signif_label controls the display text ("P=" vs "fdr="); if NULL,
    ## defaults to "P" when signif_stat=="P.Value", otherwise to
    ## signif_stat itself.
    ##
    ## Each of the three annotation pieces (base/Step1, carrier/Step3b,
    ## mediator/Step3a) is included ONLY if its columns are actually
    ## present in `stats` -- e.g. if you haven't computed Step 3b yet,
    ## just omit *_carrier columns and that piece is silently dropped
    ## rather than erroring.
    ## ------------------------------------------------------------------
    if (is.null(signif_label)) signif_label <- if (identical(signif_stat, "P.Value")) "P" else signif_stat

    stat_cols <- function(suffix_part) {
        list(pval = paste0(signif_stat, "_", anno_stat_suffix, suffix_part),
             t    = paste0("t_", anno_stat_suffix, suffix_part))
    }
    base_cols    <- stat_cols("")
    carrier_cols <- stat_cols("_carrier")
    med_cols     <- stat_cols("_med")
    mediated_col <- paste0(anno_stat_suffix, "_mediated")

    sig_stars <- function(p) dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

    stats_filter <- stats |>
        dplyr::filter(mediator == mediator_gene,
                      .data[[outcome_col]] %in% gene,
                      med_cl_test == med_clus,
                      outcome_cl == clus) |>
        dplyr::mutate(Var1 = .data[[outcome_col]])

    if (nrow(stats_filter) == 0) {
        stop(sprintf("No rows in stats for mediator '%s' matching requested outcome gene(s).", mediator_gene))
    }

    n <- nrow(stats_filter)
    has_cols <- function(cols) all(unlist(cols) %in% names(stats_filter))

    ## base/Step 1 piece -- feeds the "unadjusted" panel label
    if (has_cols(base_cols)) {
        p_ <- stats_filter[[base_cols$pval]]; t_ <- stats_filter[[base_cols$t]]
        stats_filter$anno_str_unadj <- sprintf("Base: %s=%.2e%s\nt=%.2f", signif_label, p_, sig_stars(p_), t_)
    } else {
        stats_filter$anno_str_unadj <- rep(NA_character_, n)
        message(sprintf("Note: %s/%s not found in stats -- skipping base annotation.", base_cols$pval, base_cols$t))
    }

    ## carrier/Step 3b, mediator/Step 3a, and the mediated flag -- each
    ## dropped independently if its columns aren't present
    carrier_piece <- if (has_cols(carrier_cols)) {
        p_ <- stats_filter[[carrier_cols$pval]]; t_ <- stats_filter[[carrier_cols$t]]
        sprintf("Adj for %s: %s=%.2e%s\nt=%.2f", mediator_gene, signif_label, p_, sig_stars(p_), t_)
    } else {
        message(sprintf("Note: %s/%s not found in stats -- skipping carrier-adjusted annotation.", carrier_cols$pval, carrier_cols$t))
        rep(NA_character_, n)
    }

    med_piece <- if (has_cols(med_cols)) {
        p_ <- stats_filter[[med_cols$pval]]; t_ <- stats_filter[[med_cols$t]]
        sprintf("M~Y: %s=%.2e%s\nt=%.2f", signif_label, p_, sig_stars(p_), t_)
    } else {
        message(sprintf("Note: %s/%s not found in stats -- skipping mediator-association annotation.", med_cols$pval, med_cols$t))
        rep(NA_character_, n)
    }

    mediated_piece <- if (mediated_col %in% names(stats_filter)) {
        ifelse(stats_filter[[mediated_col]] %in% TRUE, "** MEDIATED **", NA_character_)
    } else {
        message(sprintf("Note: %s not found in stats -- skipping mediated flag.", mediated_col))
        rep(NA_character_, n)
    }

    stats_filter$anno_str_adj <- vapply(seq_len(n), function(i) {
        parts <- c(carrier_piece[i], med_piece[i], mediated_piece[i])
        paste(parts[!is.na(parts)], collapse = "\n")
    }, character(1))

    ## ------------------------------------------------------------------
    ## subset both objects to their respective clusters, then to the
    ## donors present in BOTH -- same donor set feeds both panels
    ## ------------------------------------------------------------------
    sce_out <- sce[, sce[[cluster_col]] == clus]
    sce_med <- sce_mediator[, sce_mediator[[med_cluster_col]] == med_clus]

    common_brnum <- intersect(as.character(sce_out$BrNum), as.character(sce_med$BrNum))
    if (length(common_brnum) < min_donors) {
        stop(sprintf("Only %d donors shared between outcome cluster '%s' and mediator cluster '%s' (min_donors=%d).",
                     length(common_brnum), clus, med_clus, min_donors))
    }

    sce_out <- sce_out[, match(common_brnum, as.character(sce_out$BrNum))]
    sce_med <- sce_med[, match(common_brnum, as.character(sce_med$BrNum))]

    rownames(sce_med) <- rowData(sce_med)[[gene_col]]
    stopifnot(mediator_gene %in% rownames(sce_med))
    med_vec <- as.numeric(assay(sce_med[mediator_gene, ], "logcounts"))

    ## ------------------------------------------------------------------
    ## panel 1: unadjusted cleanY (same convention as plot_DEG_express)
    ## ------------------------------------------------------------------
    if (is.matrix(mod)) {
        my_mod <- mod[match(common_brnum, as.character(sce$BrNum)), ]
    } else if (inherits(mod, "formula")) {
        my_mod <- model.matrix(mod, colData(sce_out))
    }

    assays(sce_out)$cleanY <- jaffelab::cleaningY(logcounts(sce_out), mod = my_mod, P = cleanY_P)

    ## ------------------------------------------------------------------
    ## panel 2: cleanY with the mediator appended as an additional
    ## trailing covariate to regress out -- cleanY_P is unchanged since
    ## med_vec is appended AFTER the existing covariate-of-interest
    ## columns, not inserted before them.
    ## ------------------------------------------------------------------
    if (is.matrix(mod)) {
        my_mod_adj <- cbind(my_mod, med_vec = med_vec)
    } else if (inherits(mod, "formula")) {
        sce_out$med_vec <- med_vec
        my_mod_adj <- model.matrix(stats::update(mod, ~ . + med_vec), colData(sce_out))
    }

    assays(sce_out)$cleanY_adjM <- jaffelab::cleaningY(logcounts(sce_out), mod = my_mod_adj, P = cleanY_P)

    rownames(sce_out) <- rowData(sce_out)[[gene_col]]
    sce_out <- sce_out[unique(stats_filter[[outcome_col]]), ]

    ## ------------------------------------------------------------------
    ## build the two panels and combine with patchwork
    ## ------------------------------------------------------------------
    p_unadj <- plot_gene_express(
        sce = sce_out,
        genes = rownames(sce_out),
        assay_name = "cleanY",
        category = category_col,
        color_pal = color_pal,
        title = clus,
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
        ggplot2::labs(subtitle = "unadjusted") +
        ggplot2::geom_label(
            data = stats_filter, ggplot2::aes(x = -Inf, y = -Inf, label = anno_str_unadj),
            alpha = 0.5, vjust = "inward", hjust = "inward", size = 2.5
        )

    p_adj <- plot_gene_express(
        sce = sce_out,
        genes = rownames(sce_out),
        assay_name = "cleanY_adjM",
        category = category_col,
        color_pal = color_pal,
        title = clus,
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
        ggplot2::labs(subtitle = sprintf("adjusted for %s | %s", mediator_gene, med_clus)) +
        ggplot2::geom_label(
            data = stats_filter, ggplot2::aes(x = -Inf, y = -Inf, label = anno_str_adj),
            alpha = 0.5, vjust = "inward", hjust = "inward", size = 2.5
        )

    if (isTRUE(plot_mediator_panel)) {
        if (is.null(mediator_stats)) {
            stop("plot_mediator_panel=TRUE requires mediator_stats (a plot_DEG_express()-style DE table for the mediator's own cluster).")
        }
        p_med <- plot_DEG_express(
            sce = sce_mediator,
            stats = mediator_stats,
            clus = med_clus,
            gene = mediator_gene,
            pval_col = mediator_pval_col,
            fc_col = mediator_fc_col,
            gene_col = gene_col,
            cluster_col = med_cluster_col,
            category_col = category_col,
            mod = if (is.null(mediator_mod)) mod else mediator_mod,
            cleanY_P = if (is.null(mediator_cleanY_P)) cleanY_P else mediator_cleanY_P,
            color_pal = color_pal,
            plot_points = plot_points,
            ncol = ncol
        ) +
            ggplot2::labs(title = med_clus, subtitle = sprintf("mediator"))

        return(p_med + p_unadj + p_adj)
    }

    p_unadj + p_adj
}
