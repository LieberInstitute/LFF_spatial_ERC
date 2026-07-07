#' Plot outcome gene expression with and without adjusting for a mediator gene
#'
#' Adapted from `plot_DEG_express()` for the Baron & Kenny mediation
#' workflow. Where `plot_DEG_express()` shows one cleaned expression panel
#' per outcome gene (covariates regressed out), this function shows TWO
#' panels side by side for each outcome gene:
#'   - "Unadjusted": outcome gene expression with the usual nuisance
#'     covariates (Sex/Age/Anc_Afr/etc, per `mod`/`cleanY_P`) regressed out
#'     -- this is the Step 1 (X->Y) view, identical to `plot_DEG_express()`.
#'   - "Adjusted for {mediator}": the same, but with the mediator gene's
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
#' shape of `mediation_validation_details()`'s output (columns `outcome`,
#' `mediator`, `P.Value_base`, `t_base`, `P.Value_carrier`, `t_carrier`,
#' `P.Value_med`, `t_med`, `mediated`).
#' @param clus A `character(1)` outcome cluster (matches `cluster_col` in
#' `colData(sce)`).
#' @param med_clus A `character(1)` mediator cluster (matches
#' `med_cluster_col` in `colData(sce_mediator)`). This is the cluster
#' label actually used to pull mediator expression, e.g. `med_cl_test`
#' in the APOE-stratified scenarios.
#' @param mediator A `character(1)` mediator gene name (must be in
#' `rownames(sce_mediator)` after `gene_col` is applied, and must match
#' `stats$mediator`).
#' @param gene A `character()` vector of outcome gene(s) to plot (subset
#' of `stats$outcome`).
#' @param outcome_col,mediator_col The `character(1)` column names in
#' `stats` identifying the outcome gene and mediator gene respectively.
#' @param gene_col The `character(1)` name of the `rowData()` column
#' (same on both `sce` and `sce_mediator`) giving the gene symbol used to
#' match against `gene`/`mediator`.
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
#'
#' @return A `patchwork` object: two `ggplot2` panels (unadjusted,
#' adjusted-for-mediator) side by side, one facet per outcome gene within
#' each panel.
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
#'     mediator = "NPTXR",
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
                                       mediator,
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
                                       min_donors = 10) {

    stopifnot(cluster_col %in% colnames(colData(sce)))
    stopifnot(med_cluster_col %in% colnames(colData(sce_mediator)))
    stopifnot(clus %in% sce[[cluster_col]])
    stopifnot(med_clus %in% sce_mediator[[med_cluster_col]])
    stopifnot(outcome_col %in% colnames(stats))
    stopifnot(mediator_col %in% colnames(stats))
    stopifnot(all(gene %in% stats[[outcome_col]]))
    stopifnot(mediator %in% stats[[mediator_col]])

    # RCMD fix
    rank_int <- Symbol <- anno_str <- Var1 <- NULL

    ## ------------------------------------------------------------------
    ## annotation strings, one per outcome gene: base (Step 1), carrier
    ## adjusted for M (Step 3b), and M's own association (Step 3a)
    ## ------------------------------------------------------------------
    sig_stars <- function(p) dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

    stats_filter <- stats |>
        dplyr::filter(.data[[mediator_col]] == mediator, .data[[outcome_col]] %in% gene) |>
        dplyr::mutate(
            Var1 = .data[[outcome_col]],
            anno_str_unadj = sprintf("Base: P=%.2e%s\nt=%.2f", P.Value_base, sig_stars(P.Value_base), t_base),
            anno_str_adj   = sprintf("Adj for %s: P=%.2e%s\nt=%.2f\nM~Y: P=%.2e%s%s",
                                      mediator, P.Value_carrier, sig_stars(P.Value_carrier), t_carrier,
                                      P.Value_med, sig_stars(P.Value_med),
                                      ifelse(isTRUE(mediated), "\n** MEDIATED **", ""))
        )

    if (nrow(stats_filter) == 0) {
        stop(sprintf("No rows in stats for mediator '%s' matching requested outcome gene(s).", mediator))
    }

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
    stopifnot(mediator %in% rownames(sce_med))
    med_vec <- as.numeric(assay(sce_med[mediator, ], "logcounts"))

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
        title = sprintf("%s (unadjusted)", clus),
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
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
        title = sprintf("%s (adjusted for %s | %s)", clus, mediator, med_clus),
        plot_points = plot_points,
        ncol = ncol,
        plot_type = "boxplot",
        free_y = TRUE
    ) +
        ggplot2::geom_label(
            data = stats_filter, ggplot2::aes(x = -Inf, y = -Inf, label = anno_str_adj),
            alpha = 0.5, vjust = "inward", hjust = "inward", size = 2.5
        )

    p_unadj + p_adj
}
