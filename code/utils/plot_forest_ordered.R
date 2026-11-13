# Add the approximate 95% normal intervals used by crumblr::plotForest.
# These describe model effects on the CLR scale, not raw proportions.
forest_intervals <- function(x) {
    tab <- tibble::as_tibble(x)
    stopifnot(all(c("beta", "se") %in% names(tab)))
    tab$ci_lower <- tab$beta - 1.96 * tab$se
    tab$ci_upper <- tab$beta + 1.96 * tab$se
    tab
}

# cell_order is top-to-bottom. A master list can be reused for subsets:
# absent labels are ignored; unlisted plotted labels are appended in their
# original plot order. NULL preserves hierarchy order for treedata inputs
# and top-to-bottom row order for data frames/tibbles. Tabular input requires
# label, beta, and se; if n_features is present, only n_features == 1 is used.
plotForestOrdered <- function(x, cell_order = NULL, ...,
                              low = "blue", mid = "grey70", high = "red",
                              hide = FALSE) {
    if (!is.null(cell_order) &&
        (!is.character(cell_order) || anyNA(cell_order) || anyDuplicated(cell_order))) {
        stop("cell_order must be a character vector of unique, non-missing labels.")
    }
    if (is.data.frame(x)) {
        required <- c("label", "beta", "se")
        missing <- setdiff(required, names(x))
        if (length(missing)) {
            stop("Missing required columns: ", paste(missing, collapse = ", "))
        }
        tab <- tibble::as_tibble(x)
        if ("n_features" %in% names(tab)) {
            tab <- tab[!is.na(tab$n_features) & tab$n_features == 1, ]
        }
        if (!nrow(tab)) stop("No cell-type rows to plot.")
        if (anyNA(tab$label) || anyDuplicated(as.character(tab$label))) {
            stop("Cell-type labels must be unique and non-missing.")
        }
        if (!is.numeric(tab$beta) || !is.numeric(tab$se) ||
            any(!is.finite(tab$beta)) || any(!is.finite(tab$se)) || any(tab$se < 0)) {
            stop("beta and se must be finite numeric values, with se >= 0.")
        }
        tab <- forest_intervals(tab)
        tab$label <- factor(as.character(tab$label),
                            levels = rev(as.character(tab$label)))
        beta_max <- max(abs(tab$beta))
        p <- ggplot2::ggplot(tab, ggplot2::aes(label, beta)) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                color = "grey", linewidth = 1) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                                   width = 0) +
            ggplot2::geom_point(ggplot2::aes(color = beta)) +
            ggplot2::theme_classic() + ggplot2::coord_flip() +
            ggplot2::labs(x = "", y = "Effect size") +
            ggplot2::scale_color_gradient2(
                name = expression(beta), low = low, mid = mid, high = high,
                midpoint = 0, limits = c(-beta_max, beta_max))
        if (hide) {
            p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                    axis.ticks.y = ggplot2::element_blank(),
                                    legend.position = "bottom")
        }
    } else {
        p <- crumblr::plotForest(x, ..., low = low, mid = mid, high = high,
                                hide = hide)
        p$data <- forest_intervals(p$data)
    }
    if (!is.null(cell_order)) {
        plotted <- as.character(p$data$label)
        original_order <- rev(levels(p$data$label))
        ordered <- c(intersect(cell_order, plotted),
                     setdiff(original_order, cell_order))
        # crumblr uses coord_flip(), so its discrete x axis is vertical.
        p <- p + ggplot2::scale_x_discrete(limits = rev(ordered))
    }
    p
}

# Export single-coefficient results with moderated SEs and model-based 95% CIs.
# Match by row name because topTable sorts results by significance.
topTableWithSE <- function(fit, coef) {
    stopifnot(length(coef) == 1L)
    tab <- variancePartition::topTable(fit, coef = coef, number = Inf,
                                     confint = 0.95)
    idx <- match(rownames(tab), rownames(fit$coefficients))
    stopifnot(!anyNA(idx))
    posterior_variance <- fit$s2.post
    if (is.matrix(posterior_variance)) {
        posterior_variance <- posterior_variance[, coef]
    }
    tab$se <- fit$stdev.unscaled[idx, coef] * sqrt(posterior_variance[idx])
    # Also supply the normal-approximation endpoints used by plotForest.
    tab$forest_ci_lower <- tab$logFC - 1.96 * tab$se
    tab$forest_ci_upper <- tab$logFC + 1.96 * tab$se
    tab
}
