library(here)
library(tidyverse)
library(sessioninfo)

in_dir = here("processed-data", "14_MOFA", "08_project_xenium")
plot_dir = here("plots", "14_MOFA", "09_plot_xenium")
project_colors_path = here("processed-data", "project_colors.Rdata")
tau_colors = c(`t-` = "#684F7D", `t+` = "#AFA4B6")

dir.create(plot_dir, showWarnings = FALSE)
load(project_colors_path, verbose = TRUE)

################################################################################
#   Functions
################################################################################

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
        labs(x = "Samples", y = "Factor 3 Weight") +
        theme_bw(base_size = 12) +
        theme(legend.position = "None")  + 
        geom_label(
            data = t_df |>
                mutate(
                    p_value_label = sprintf(
                        "p=%.2e%s", p_value, ifelse(p_value < 0.05, "*", "")
                    )
                ),
            aes(x = Inf, y = Inf, label = p_value_label),
            size = 4, vjust = "inward", hjust = "inward"
        )
    pdf(file.path(plot_dir, sprintf("%s_boxplot.pdf", covariate)), width = 10)
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

covariate_boxplot(
    factor_df = factor_df,
    t_df = factor_t_test(factor_df, covariate = "APOE_carrier"),
    covariate = "APOE_carrier",
    covariate_colors = APOE_carrier_colors
)

covariate_boxplot(
    factor_df = factor_df,
    t_df = factor_t_test(factor_df, covariate = "taupathy"),
    covariate = "taupathy",
    covariate_colors = tau_colors
)

session_info()
