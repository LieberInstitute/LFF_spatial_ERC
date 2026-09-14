library(here)
library(tidyverse)
library(sessioninfo)

in_dir = here("processed-data", "14_MOFA", "08_project_xenium")
plot_dir = here("plots", "14_MOFA", "09_plot_xenium")
out_path = here(
    "processed-data", "14_MOFA", "09_plot_xenium", "factor_t_test_results.csv"
)
project_colors_path = here("processed-data", "project_colors.Rdata")
tau_colors = c(`t-` = "#684F7D", `t+` = "#AFA4B6")

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)
load(project_colors_path, verbose = TRUE)

################################################################################
#   Functions
################################################################################

#   Format a (typically FDR-adjusted) p-value for plot labels: scientific
#   notation below 0.001, otherwise 3 decimal places - avoids round() collapsing
#   small but meaningfully-different p-values down to "0" or "0.00"
format_fdr = function(x) {
    ifelse(x < 0.001, sprintf("%.1e", x), sprintf("%.3f", x))
}

factor_t_test = function(factor_df, covariate) {
    t_df_list = list()
    for (this_facet_label in unique(factor_df$facet_label)) {
        t_result = t.test(
            as.formula(sprintf("xen_factor_score ~ %s", covariate)),
            data = factor_df |>
                filter(facet_label == this_facet_label)
        )
        t_df_list[[this_facet_label]] = tibble(
            facet_label = this_facet_label,
            t_statistic = t_result$statistic,
            p_value = t_result$p.value
        )
    }
    t_df = bind_rows(t_df_list)

    return(t_df)
}

covariate_boxplot = function(factor_df, t_df, covariate, covariate_colors) {
    p = ggplot(factor_df) +
        geom_boxplot(
            aes(x = !!sym(covariate), y = xen_factor_score, fill = !!sym(covariate)),
            outlier.shape = NA
        ) +
        geom_jitter(
            aes(x = !!sym(covariate), y = xen_factor_score, fill = !!sym(covariate)),
            width = 0.2, size = 2
        ) +
        facet_wrap(~facet_label, nrow = 3, scales = "free_y") +
        scale_fill_manual(values = covariate_colors) +
        labs(x = "Samples", y = "Factor 4 Weight") +
        theme_bw(base_size = 12) +
        theme(legend.position = "None")  + 
        geom_label(
            data = t_df |>
                mutate(
                    p_value_label = sprintf(
                        "FDR=%s%s", format_fdr(fdr), ifelse(fdr < 0.05, "*", "")
                    )
                ),
            aes(x = Inf, y = Inf, label = p_value_label),
            size = 4, vjust = "inward", hjust = "inward"
        )
    pdf(file.path(plot_dir, sprintf("%s_boxplot.pdf", covariate)), width = 11)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

factor_df = list.files(
        in_dir, pattern = "factor.*\\.csv$", full.names = TRUE
    ) |>
    map_dfr(read_csv, show_col_types = FALSE) |>
    mutate(facet_label = sprintf('SpD: %s; Astro: %s', domain, astro_group))

t_df = rbind(
    factor_t_test(factor_df, covariate = "APOE_carrier") |>
        mutate(covariate = "APOE_carrier"),
    factor_t_test(factor_df, covariate = "taupathy") |>
        mutate(covariate = "taupathy")
)

#   FDR-correct within each covariate's family of tests (across all SpD/astro
#   facets tested for that covariate) rather than pooling both covariates
t_df = t_df |>
    group_by(covariate) |>
    mutate(fdr = p.adjust(p_value, method = "BH")) |>
    ungroup()

covariate_boxplot(
    factor_df = factor_df,
    t_df = t_df |> filter(covariate == "APOE_carrier"),
    covariate = "APOE_carrier",
    covariate_colors = APOE_carrier_colors
)

covariate_boxplot(
    factor_df = factor_df,
    t_df = t_df |> filter(covariate == "taupathy"),
    covariate = "taupathy",
    covariate_colors = tau_colors
)

#   Also save t-stat, raw p-value, and FDR-adjusted p-value for each covariate
#   for later (FDR computed within each covariate's family of tests above)
t_df |>
    mutate(
        SpD_subset = str_extract(
            facet_label, "^SpD: ([^;]+);", group = 1
        ),
        astro_subset = str_extract(
            facet_label, "Astro: (.*)$", group = 1
        )
    ) |>
    dplyr::rename(t_stat = t_statistic, p_val = p_value) |>
    write_csv(out_path)

session_info()
