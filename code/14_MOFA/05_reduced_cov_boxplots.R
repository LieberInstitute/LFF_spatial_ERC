library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)

model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
plot_path = here("plots", "14_MOFA", "05_reduced_cov_boxplots", "boxplot.pdf")
project_colors_path = here("processed-data", "project_colors.Rdata")
specific_factor = "Factor3"
tau_colors = c(`t-` = "#684F7D", `t+` = "#AFA4B6")

dir.create(dirname(plot_path), showWarnings = FALSE)

model = load_model(model_path)
load(project_colors_path, verbose = TRUE)

factor_df = get_tidy_factors(
    model = model,
    metadata = samples_metadata(model),
    factor = specific_factor,
    sample_id_column = "sample"
)

test_vars = c('APOE_carrier', "taupathy")
names(test_vars) = test_vars

assoc_list = map(
    test_vars,
    ~ get_associations(
        model = model,
        metadata = samples_metadata(model),
        sample_id_column = "sample",
        test_variable = .x,
        test_type = "categorical",
        group = FALSE
    )
)

assoc_tb = do.call("bind_rows", assoc_list) |>
    mutate(
        signif = case_when(
            adj_pvalue < 0.001 ~ "***",
            adj_pvalue < 0.01 ~ "**",
            adj_pvalue < 0.05 ~ "*",
            TRUE ~ ""
        )
    )

p = factor_df |>
    dplyr::rename(factor_val = value) |>
    select(sample, APOE_carrier, taupathy, factor_val) |>
    pivot_longer(!c(sample, factor_val), names_to = "term") |>
    mutate(term = factor(term, levels = test_vars)) |>
    ggplot() +
        geom_boxplot(
            aes(x = value, y = factor_val, fill = value), outlier.shape = NA
        ) +
        geom_jitter(
            aes(x = value, y = factor_val, fill = value), width = 0.2, size = 3
        ) +
        facet_wrap(~term, nrow = 1, scales = "free_x") +
        scale_fill_manual(values = c(APOE_carrier_colors, tau_colors)) +
        labs(x = "Samples", y = "Factor 3 Weight") +
        theme_bw(base_size = 30) +
        theme(legend.position = "None")  + 
        geom_label(
            data = assoc_tb |> 
                filter(Factor == "Factor3", term %in% test_vars) |>
                mutate(term = factor(term, levels = test_vars)),
            aes(x = Inf, y = Inf, label = sprintf("FDR=%.2e%s", adj_pvalue, signif)),
            size = 8,
            vjust = "inward", 
            hjust = "inward"
        )
pdf(plot_path, height = 5)
print(p)
dev.off()

session_info()
