## Draft: scenario-specific Xenium evidence for the best-supported pairs.
## Run from the project root. Does not change the existing Sankey script.
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(here)
})

out_dir <- here("plots", "22_Mediation", "05_Mediation_general_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results <- read.csv(here("processed-data", "22_Mediation", "03_Mediation_Xenium",
                         "mediation_summary-all_scenarios_Pval0.10.csv"))
threshold <- 0.10

# Select pairs, then retain ALL their tested scenarios.
selected <- results |>
    filter(base_valid, mediated) |>
    distinct(med_cl, mediator, outcome)
evidence <- results |>
    semi_join(selected, by = c("med_cl", "mediator", "outcome"))

# Load mediator DE for every scenario, including those absent from the hits CSV.
mediator_stats <- bind_rows(lapply(unique(evidence$mediator_datatype), function(datatype) {
    de <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE",
                      datatype, paste0("DGE_results_carrier_", datatype, ".Rds")))
    if (grepl("_SpX$", datatype)) de$cluster <- de$cluster_xSpD
    de |>
        transmute(mediator_datatype = datatype, med_cl_test = cluster,
                  mediator = gene_name, mediator_p = vlmf_P.Value)
}))
keys <- c("mediator_datatype", "med_cl_test", "mediator")
stopifnot(!anyDuplicated(mediator_stats[keys]))

support_levels <- c("All screening criteria", "Conditional association; no mediator DE",
                    "Loss of significance only", "Other tested result")
palette <- setNames(c("#147D83", "#4988BD", "#DCA444", "#CCD1D5"), support_levels)

evidence <- evidence |>
    left_join(mediator_stats, by = keys) |>
    mutate(
        mediator_de = !is.na(mediator_p) & mediator_p < threshold,
        support = case_when(
            base_valid & mediated & mediator_de ~ support_levels[1],
            base_valid & mediated ~ support_levels[2],
            base_valid & !med_sig ~ support_levels[3],
            TRUE ~ support_levels[4]
        ),
        scenario = case_when(
            grepl("_SpX$", mediator_datatype) ~ sub("^.*_x", "", med_cl_test),
            grepl("APOE_high", med_cl_test) ~ "APOE high",
            grepl("APOE_low", med_cl_test) ~ "APOE low",
            grepl("^nn_", outcome_cl) ~ "Astro-adjacent",
            outcome_cl == "Oligo.3" ~ "Overall",
            TRUE ~ NA_character_
        ),
        pair = paste0(med_cl, "  ", mediator, "  →  ", outcome),
        conditional_strength = -log10(pmax(P.Value_med_vec, .Machine$double.xmin)),
        support = factor(support, levels = support_levels)
    )
stopifnot(!anyNA(evidence$scenario), !anyNA(evidence$mediator_p),
          !anyDuplicated(evidence[c("pair", "scenario")]))

scenario_key <- tibble(
    scenario = c("Overall", "Astro-adjacent", "GL", "L1_inhib", "L2", "L5", "L6",
                 "LD_glia", "Vasc_adj", "WMd", "WMuf", "APOE high", "APOE low"),
    family = c(rep("Population", 2), rep("Spatial domain", 9), rep("APOE neighborhood", 2))
)
pair_order <- evidence |>
    group_by(pair) |>
    summarise(rank = min(as.integer(support)), .groups = "drop") |>
    arrange(rank, pair) |>
    pull(pair)

# Blank cells represent no result in the summary, not failed validation.
grid <- expand_grid(pair = pair_order, scenario = scenario_key$scenario) |>
    left_join(scenario_key, by = "scenario") |>
    mutate(pair = factor(pair, levels = rev(pair_order)),
           scenario = factor(scenario, levels = scenario_key$scenario),
           family = factor(family, levels = unique(scenario_key$family)))
evidence <- evidence |>
    left_join(scenario_key, by = "scenario") |>
    mutate(pair = factor(pair, levels = rev(pair_order)),
           scenario = factor(scenario, levels = scenario_key$scenario),
           family = factor(family, levels = unique(scenario_key$family)))

p <- ggplot(grid, aes(scenario, pair)) +
    geom_tile(fill = "#F6F7F8", color = "white", linewidth = 1) +
    geom_point(data = evidence,
               aes(size = conditional_strength, fill = support),
               shape = 21, color = "white", stroke = 0.6) +
    # A small plus marks mediator DE independently of the support category.
    geom_point(data = filter(evidence, mediator_de),
               shape = 3, size = 1.6, stroke = 0.5, color = "#172C35") +
    facet_grid(~family, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_size_continuous(range = c(2, 10), limits = c(0, NA),
                          breaks = c(0, 1, 2, 3)) +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x), drop = TRUE) +
    labs(title = "Xenium support for candidate mediation pairs",
         subtitle = "Pairs retained if baseline replication and conditional mediator association coincide with a nonsignificant adjusted APOE effect",
         x = NULL, y = NULL,
         fill = "Within-scenario evidence", size = expression(-log[10](P[M~"|"~APOE])),
         caption = paste("+  Mediator DE P < .10.  Dot size: conditional mediator–outcome association.  Blank: no result available.",
                         "All screening thresholds use nominal P < .10. Scenarios reuse donors and are not independent replications.",
                         "Loss of significance does not itself establish attenuation or a causal indirect effect.", sep = "\n")) +
    guides(fill = guide_legend(order = 1, nrow = 2, byrow = TRUE,
                               override.aes = list(size = 4)),
           size = guide_legend(order = 2, nrow = 1,
                               override.aes = list(fill = "#778A96", color = "#778A96"))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0),
          axis.text.y = element_text(color = "#172C35", face = "bold"),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
          plot.title = element_text(size = 19, face = "bold"),
          plot.subtitle = element_text(size = 10, margin = margin(b = 18)),
          plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.3),
          legend.position = "bottom", legend.box = "vertical",
          plot.margin = margin(18, 18, 12, 18))

ggsave(file.path(out_dir, "mediation_xenium_support_dotplot_draft.png"), p,
       width = 14, height = 6.8, dpi = 220, bg = "white")
write.csv(evidence, file.path(out_dir, "mediation_xenium_support_dotplot_data.csv"), row.names = FALSE)
print(evidence |> count(pair, support))
