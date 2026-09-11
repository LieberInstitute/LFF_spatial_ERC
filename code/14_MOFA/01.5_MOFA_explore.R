## Louise Huuki-Myers, July 2025
## Multicellular factor analysis - exploratory / covariate-specific plots
## Split off from 01_MOFA.R: these plots are more likely to change between
## iterations, so this script loads the already-fit MOFA model (rather than
## re-running run_mofa()) and re-derives what it needs from 01_MOFA.R's saved
## outputs. Run 01_MOFA.R first for a given --datatype before running this.

#### set up ####
library("SpatialExperiment")
library("MOFAcellulaR")
library("MOFA2")
library("tidyverse")
library("GGally")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"

data_dir <- here("processed-data", "14_MOFA", "01_MOFA", opt$datatype)
plot_dir <- here("plots", "14_MOFA", "01_MOFA", opt$datatype)

## colors
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)
load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE)
load(here("processed-data", "SpD_colors.Rdata"), verbose = TRUE)

#### Load fitted MOFA model + outputs from 01_MOFA.R ####
message(Sys.time(), " - Load MOFA model")
out_path <- here(data_dir, "model.hdf5")
model <- load_model(out_path)

factor_df <- read_csv(here(data_dir, sprintf("MOFA_factor_df-%s.csv", opt$datatype)))
assoc_tb <- read_csv(here(data_dir, sprintf("MOFA_factor_associations_df-%s.csv", opt$datatype))) |>
    mutate(signif = case_when(adj_pvalue < 0.001 ~ "***",
                              adj_pvalue < 0.01 ~ "**",
                              adj_pvalue < 0.05 ~ "*",
                              TRUE ~ "")
    )

## assoc_list is a named list (not a flat table), so it isn't saved to csv by
## 01_MOFA.R - cheaply re-derive it here from the loaded model
test_vars <- c('APOE_carrier', 'Ancestry', "taupathy", "Sex", "Age", "Braak", "CERAD")
names(test_vars) <- test_vars

assoc_list <- map(test_vars, ~get_associations(model = model,
                                               metadata = samples_metadata(model),
                                               sample_id_column = "sample",
                                               test_variable = .x,
                                               test_type = "categorical",
                                               group = FALSE))

## gene id/name mapping saved by 01_MOFA.R (from rowData(spe))
# rd <- read_rds(here(data_dir, sprintf("MOFA_gene_rowData_%s.rds", opt$datatype)))

## factor 4 boxplots

# tau_colors <- c(`t-` = "#493657", `t+` = "#BF2626")
tau_colors <- c(`t-` = "#684F7D", `t+` = "#AFA4B6")

F4_weights_boxplot <- factor_df |>
    filter(Factor == "Factor4") |>
    select(sample, APOE_carrier, taupathy, Sex, Factor4 = value) |>
    pivot_longer(!c(sample, Factor4), names_to = "term") |>
    mutate(term = factor(term, levels = c("APOE_carrier", "taupathy", "Sex")))|>
    ggplot() +
    geom_boxplot(aes(x = value, y = Factor4, fill = value), outlier.shape = NA) +
    geom_jitter(aes(x = value, y = Factor4, fill = value), width = 0.1) +
    facet_wrap(~term, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(APOE_carrier_colors, sex_colors, tau_colors)) +
    labs(x = "samples", y = "factor 4 Weight") +
    theme_bw() +
    theme(legend.position = "None")  + 
    ggplot2::geom_label(
        data = assoc_tb |> 
            filter(Factor == "Factor4", 
                   term %in% c("APOE_carrier", "taupathy", "Sex")) |>
            mutate(term = factor(term, levels = c("APOE_carrier", "taupathy", "Sex"))), 
        ggplot2::aes(x = Inf, y = Inf, label = sprintf("FDR=%.2e%s", adj_pvalue, signif)),
        size = 2,
        vjust = "inward", 
        hjust = "inward"
    )

ggsave(F4_weights_boxplot, filename = here(plot_dir, "Factor4_weights_boxplot.png"), height = 4, width = 4)

# NOTE (flagged, not removed): this block referenced `weights_boxplot` and `var`,
# which only exist as local variables inside the factor_boxplot() function
# defined in 01_MOFA.R - they were never in scope here, so this would error
# ("object 'weights_boxplot' not found") in the original combined script too.
# Commented out pending a decision on what this was meant to do.
# weights_boxplot <- weights_boxplot + ggplot2::geom_label(
#     data = assoc_tb |> filter(term == var),
#     ggplot2::aes(x = -Inf, y = -Inf, label = sprintf("pval=%.2e%s", adj_pvalue, signif)),
#     alpha = 0.5,
#     vjust = "inward",
#     hjust = "inward",
#     size = 2.5
# )


carrier_tau_colors <- c(`E2+ t-` = "#398A84",
                        `E2+ t+` = "#7FBFBB",
                        `E4+ t-` = "#D46B43",
                        `E4+ t+` = "#DF9E85")

F4_carrier_tau_boxplot <- factor_df |>
    filter(Factor == "Factor4") |>
    mutate(carrier_tau = paste(APOE_carrier, taupathy)) |>
    ggplot() +
    geom_boxplot(aes(x = carrier_tau, y = value, fill = carrier_tau), outlier.shape = NA) +
    geom_jitter(aes(x = carrier_tau, y = value, fill = carrier_tau, shape = Sex), width = 0.1) +
    scale_fill_manual(values = carrier_tau_colors) +
    labs(x="APOE Carrier + Tau", y = "Factor4 Weight") +
    guides(fill = "none") +
    theme_bw() +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(F4_carrier_tau_boxplot, filename = here(plot_dir, "Factor4_weights_boxplot_carrier_tau_boxplot.png"), height = 4, width = 6)

ggsave(F4_carrier_tau_boxplot, filename = here(plot_dir, "Factor4_weights_boxplot_carrier_tau_boxplot_small.png"), height = 4, width = 3)

## factor 4 vs. Age

F4_df <- factor_df |>
    filter(Factor == "Factor4") |>
    mutate(carrier_tau = paste(APOE_carrier, taupathy))

F4_age_scatter <- F4_df  |>
    ggplot(aes(x = Age, y = value, color = carrier_tau, shape = Sex)) +
    scale_color_manual(values = carrier_tau_colors) +
    geom_point() +
    theme_bw()  +
    labs(x="Age", y = "Factor4 Weight")

ggsave(F4_age_scatter, filename = here(plot_dir, "Factor4_weights_age_scatter.png"), height = 4, width = 6)
ggsave(F4_age_scatter + facet_grid(Ancestry~APOE_carrier) , filename = here(plot_dir, "Factor4_weights_age_scatter_facet.png"), height = 4, width = 6)
ggsave(F4_age_scatter + facet_grid(Sex~APOE_carrier) , filename = here(plot_dir, "Factor4_weights_age_scatter_facet_sex.png"), height = 4, width = 6)
ggsave(F4_age_scatter + ggrepel::geom_text_repel(aes(label = sample)), filename = here(plot_dir, "Factor4_weights_age_scatter_text.pdf"), height = 4, width = 6)


F4_cor_stats <- F4_df |>
    group_by(APOE_carrier) |>
    summarise(broom::tidy(cor.test(Age, value)), .groups = "drop") |>
    dplyr::select(APOE_carrier, estimate, p.value)

# APOE_carrier estimate p.value
# <chr>           <dbl>   <dbl>
# 1 E2+             0.644  0.0130
# 2 E4+             0.484  0.0574

F4_age_scatter_fit <- F4_df |>
    ggplot(aes(x = Age, y = value, color = APOE_carrier, shape = taupathy)) +
    geom_point() +
    geom_smooth(method = "lm", aes(group = APOE_carrier, fill = APOE_carrier))+
    ggpubr::stat_cor(
        aes(group = APOE_carrier, color = APOE_carrier),
        method = "pearson",
        label.x.npc = "left",
        label.y.npc = 0.25,  
        show.legend = FALSE
    ) +
    ggpubr::stat_regline_equation(
        aes(group = APOE_carrier, color = APOE_carrier),
        label.x.npc = "left",
        label.y.npc = 0.1,  
        show.legend = FALSE
    ) +
    scale_color_manual(values = APOE_carrier_colors) +
    scale_fill_manual(values = APOE_carrier_colors) +
    theme_bw()  +
    labs(x="Age", y = "Factor4 Weight")

ggsave(F4_age_scatter_fit, filename = here(plot_dir, "Factor4_weights_age_scatter_fit.png"), height = 5, width = 6)


####  Age x taupathy interaction
F4_taupathy_cor <- F4_df |>
    group_by(taupathy) |>
    summarise(broom::tidy(cor.test(Age, value)), .groups = "drop") |>
    dplyr::select(taupathy, estimate, p.value)

F4_taupathy_interaction <- lm(value ~ Age * taupathy, data = F4_df)
broom::tidy(F4_taupathy_interaction)

# term           estimate std.error statistic p.value
# <chr>             <dbl>     <dbl>     <dbl>   <dbl>
# 1 (Intercept)    -0.274     0.348      -0.789   0.437
# 2 Age             0.00329   0.00682     0.482   0.634
# 3 taupathyt+      0.0941    0.773       0.122   0.904
# 4 Age:taupathyt+  0.00456   0.0138      0.330   0.744

F4_age_scatter_fit_tau <- factor_df |>
    filter(Factor == "Factor4") |>
    ggplot(aes(x = Age, y = value, color = taupathy, shape = APOE_carrier)) +
    geom_point() +
    geom_smooth(method = "lm", aes(group = taupathy, fill = taupathy))+
    ggpubr::stat_cor(
        aes(group = taupathy),
        method = "pearson",
        label.x.npc = "left",
        label.y.npc = 0.3, 
        show.legend = FALSE
    )  +
    ggpubr::stat_regline_equation(
        aes(group = taupathy, color = taupathy),
        label.x.npc = "left",
        label.y.npc = 0.2,  
        show.legend = FALSE
    ) +
    scale_color_manual(values = tau_colors) +
    scale_fill_manual(values = tau_colors) +
    theme_bw()  +
    labs(x="Age", y = "Factor4 Weight")

ggsave(F4_age_scatter_fit_tau, filename = here(plot_dir, "Factor4_weights_age_scatter_fit_tau.png"), height = 5, width = 6)

# slurmjobs::job_single('01.5_MOFA_explore_fine', create_shell = TRUE, memory = '10G', command = "Rscript 01.5_MOFA_explore.R --datatype sn_fine")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()