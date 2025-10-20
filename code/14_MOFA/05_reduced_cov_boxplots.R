library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

model_path = here(
    "processed-data", "14_MOFA", "01_MOFA", "sn_fine", "model.hdf5"
)
plot_path = here("plots", "14_MOFA", "05_reduced_cov_boxplots", "boxplot.pdf")
specific_factor = "Factor3"

dir.create(dirname(plot_path), showWarnings = FALSE)

model = load_model(model_path)

factor_df = get_tidy_factors(
    model = model,
    metadata = samples_metadata(model),
    factor = "all",
    sample_id_column = "sample"
)

test_vars <- c('APOE_carrier', 'Ancestry', "taupathy", "Sex", "Age","Braak","CERAD")
names(test_vars) <- test_vars

assoc_list <- map(test_vars, ~get_associations(model = model,
                                               metadata = samples_metadata(model),
                                               sample_id_column = "sample",
                                               test_variable = .x,
                                               test_type = "categorical",
                                               group = FALSE))

assoc_tb <- do.call("bind_rows", assoc_list) |>
    mutate(signif = case_when(adj_pvalue < 0.001 ~ "***",
                              adj_pvalue < 0.01 ~ "**",
                              adj_pvalue < 0.05 ~ "*",
                              TRUE ~ "")
    ) 

tau_colors <- c(`t-` = "#684F7D", `t+` = "#AFA4B6")

F3_weights_boxplot <- factor_df |>
    filter(Factor == "Factor3") |>
    select(sample, APOE_carrier, taupathy, Sex, Factor3 = value) |>
    pivot_longer(!c(sample, Factor3), names_to = "term") |>
    mutate(term = factor(term, levels = c("APOE_carrier", "taupathy", "Sex")))|>
    ggplot() +
    geom_boxplot(aes(x = value, y = Factor3, fill = value), outlier.shape = NA) +
    geom_jitter(aes(x = value, y = Factor3, fill = value), width = 0.1) +
    facet_wrap(~term, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(APOE_carrier_colors, sex_colors, tau_colors)) +
    labs(x = "samples", y = "Factor 3 Weight") +
    theme_bw() +
    theme(legend.position = "None")  + 
    ggplot2::geom_label(
        data = assoc_tb |> 
            filter(Factor == "Factor3", 
                   term %in% c("APOE_carrier", "taupathy", "Sex")) |>
            mutate(term = factor(term, levels = c("APOE_carrier", "taupathy", "Sex"))), 
        ggplot2::aes(x = Inf, y = Inf, label = sprintf("FDR=%.2e%s", adj_pvalue, signif)),
        size = 2,
        vjust = "inward", 
        hjust = "inward"
    )
