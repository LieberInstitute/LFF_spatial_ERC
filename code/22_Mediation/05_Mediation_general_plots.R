## Louise Huuki-Myers, July 2026
## plot scdotplots, snakey, and upset plots for mediation genes 


#### Set up ####

library("tidyverse")
library("here")
library("ggalluvial")
library("scDotPlot")
library("SingleCellExperiment")
library("UpSetR")
library("ComplexUpset")

plot_dir <- here("plots", "22_Mediation", "05_Mediation_general_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "00_project_prep", "cell_type_colors.V2.Rdata"), verbose = TRUE) 

astro_colors <- cell_type_colors$anno[grep("Astro", names(cell_type_colors$anno))]

#### Load Data ####
mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(pair = paste0(mediator, "|", outcome))


mediator_outcome |> dplyr::count(mediator)

mediator_outcome |> 
    dplyr::count(med_cl, mediator) |>
    arrange(-n)

# med_cl mediator   n
# 1 Astro.3    CHRM3 205
# 2 Astro.2     SV2B 199
# 3 Astro.2     FZD8 125
# 4 Astro.2    NETO1  61
# 5 Astro.1    NPTXR  22
# 6 Astro.3    IL6ST   4
# 7 Astro.1   PLPPR4   1
# 8 Astro.2    ABCA8   1
# 9 Astro.2     ST18   1

outcome_summary <- mediator_outcome |> 
    group_by(med_cl, outcome) |>
    summarise(n_med = n(),
              mediators = paste0(mediator, collapse = ", "))


outcome_summary |> filter(mediators == "CHRM3") |> dplyr::count(med_cl)
outcome_summary |> dplyr::count(mediators) ## 14 sets

outcome_summary |> filter(mediators == "FZD8, NETO1, SV2B") |> pull(outcome)
outcome_summary |> filter(grepl("SOX", outcome))
# med_cl  outcome n_med mediators        
# <chr>   <chr>   <int> <chr>            
# 1 Astro.1 SOX6        1 NPTXR            
# 2 Astro.2 SOX5        3 FZD8, NETO1, SV2B
# 3 Astro.3 SOX5        1 CHRM3  
    
outcome_summary |> filter(outcome == "CABLES1")

#### Upset plots ####
mediator_outcome_wide <- mediator_outcome |> 
    mutate(mediator = paste0(med_cl, "_", mediator)) |>
    select(mediator, outcome) |>
    mutate(mediated = TRUE) |>
    pivot_wider(names_from = "mediator", values_from = "mediated") |>
    mutate(
        across(everything(), ~replace_na(.x, FALSE))
    )

all_mediator_cl <- mediator_outcome |> 
    mutate(mediator = paste0(med_cl, "_", mediator)) |>
    distinct(mediator) |>
    pull(mediator)

mediator_outcome_wide |> filter(outcome == "SOX5")

mediator_outcome_wide |> filter(Astro.2_ABCA8)

astro_set_colors <- map(all_mediator_cl, ~upset_query(set=c(.x), fill=astro_colors[[jaffelab::ss(.x, "_")]]))

sox_queries <- list(upset_query(intersect = c("Astro.3_CHRM3", "Astro.2_FZD8", "Astro.2_NETO1", "Astro.2_SV2B"), fill="red"), #SOX5
                    upset_query(intersect = c("Astro.3_CHRM3", "Astro.2_FZD8", "Astro.2_SV2B"), fill="darkred"), #CABLES1
                    upset_query(intersect = c("Astro.1_NPTXR"), fill="seagreen")) #SOX6

my_queries <- c(astro_set_colors, sox_queries)

mediation_upsets <- ComplexUpset::upset(data = mediator_outcome_wide, 
                                        intersect = all_mediator_cl,
                                        name = "mediator",
                                        queries=my_queries,
                                        width_ratio=0.2) 

ggsave(mediation_upsets,filename =  here(plot_dir, "Mediation_gene_set_upset.pdf"), height = 5, width = 6)
ggsave(mediation_upsets,filename =  here(plot_dir, "Mediation_gene_set_upset.png"), height = 5, width = 6)


#### Snakey plots ####

#### Xenium validation support ####

support_levels <- c(
    "All screening criteria supported",
    "Conditional association; mediator not DE",
    "Loss of significance only",
    "Other tested result"
)

support_colors <- setNames(
    c("#167D8D", "#619CDB", "#E5AD45", "#E5E5E5"),
    support_levels
)

# Full results preserve tested scenarios that did not enter the hits file.
xenium_all <- read_csv(here(
    "processed-data", "22_Mediation", "03_Mediation_Xenium",
    "mediation_summary-all_scenarios_Pval0.10.csv"
))

xenium_hits <- read_csv(here(
    "processed-data", "22_Mediation", "03_Mediation_Xenium",
    "Xenium_mediation_hits.csv"
))

# Mediator DE depends on its datatype, cluster, and gene—not outcome.
# The hits file supplies these statistics for all candidate hit scenarios.
mediator_stats <- xenium_hits |>
    distinct(
        mediator_datatype, med_cl_test, mediator,
        mediatorDE_P.Value
    )

stopifnot(!anyDuplicated(
    mediator_stats[c("mediator_datatype", "med_cl_test", "mediator")]
))

xenium_support <- xenium_all |>
    left_join(
        mediator_stats,
        by = c("mediator_datatype", "med_cl_test", "mediator")
    ) |>
    mutate(
        baseline_supported = coalesce(base_valid, FALSE),
        adjusted_nonsig = !is.na(P.Value_med) & P.Value_med >= 0.10,
        conditional_supported =
            !is.na(P.Value_med_vec) & P.Value_med_vec < 0.10,
        mediator_supported =
            !is.na(mediatorDE_P.Value) & mediatorDE_P.Value < 0.10,
        
        # These flags must coincide in the SAME row/scenario.
        support = case_when(
            baseline_supported & adjusted_nonsig &
                conditional_supported & mediator_supported ~
                support_levels[1],
            
            baseline_supported & adjusted_nonsig &
                conditional_supported ~
                support_levels[2],
            
            baseline_supported & adjusted_nonsig ~
                support_levels[3],
            
            TRUE ~ support_levels[4]
        ),
        support = factor(support, levels = support_levels),
        pair_label = paste0(med_cl, ": ", mediator, " → ", outcome),
        
        # Keep full identities to avoid combining distinct analyses.
        scenario = paste(
            mediator_datatype, med_cl_test,
            outcome_datatype, outcome_cl,
            sep = " | "
        ),
        scenario_label = paste0(med_cl_test, " → ", outcome_cl)
    )

# Catch missing mediator statistics where they affect classification.
stopifnot(!any(
    xenium_support$baseline_supported &
        xenium_support$adjusted_nonsig &
        is.na(xenium_support$mediatorDE_P.Value)
))

# Show pairs with at least one candidate-support scenario.
# Remove this filter to display all 58 tested pairs.
candidate_pairs <- xenium_support |>
    filter(support != support_levels[4]) |>
    distinct(med_cl, mediator, outcome)

candidate_support <- xenium_support |>
    semi_join(candidate_pairs, by = c("med_cl", "mediator", "outcome"))

#### Sankey: one equal-width link per unique pair ####

pair_support <- candidate_support |>
    group_by(med_cl, mediator, outcome, pair_label) |>
    summarise(
        best_rank = min(as.integer(support)),
        n_tested = n_distinct(scenario),
        n_all_criteria = sum(support == support_levels[1]),
        .groups = "drop"
    ) |>
    mutate(
        support = factor(
            support_levels[best_rank],
            levels = support_levels
        )
    )

mediation_xenium_support_sankey <- ggplot(
    pair_support,
    aes(axis1 = mediator, axis2 = outcome, y = 1)
) +
    geom_alluvium(aes(fill = support), width = 1 / 4, alpha = 0.85) +
    geom_stratum(width = 1 / 4, fill = "white", color = "grey50") +
    geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 3
    ) +
    scale_fill_manual(values = support_colors, drop = FALSE) +
    scale_x_discrete(
        limits = c("Mediator", "Outcome"),
        expand = c(0.10, 0.10)
    ) +
    facet_wrap(~med_cl, scales = "free_y", nrow = 2) +
    labs(
        fill = "Strongest support\nwithin one scenario"
    ) +
    theme_void() +
    theme(
        legend.position = "bottom",
        plot.caption = element_text(size = 9)
    ) +
    guides(fill = guide_legend(ncol = 1))

ggsave(
    here(plot_dir, "mediation_xenium_support_sankey.png"),
    mediation_xenium_support_sankey,
    width = 4, height = 4, dpi = 300
)

#### Evidence matrix: retain every tested scenario for these pairs ####

pair_order <- pair_support |>
    arrange(best_rank, med_cl, mediator, outcome) |>
    pull(pair_label)

scenario_key <- candidate_support |>
    distinct(scenario, scenario_label) |>
    arrange(scenario)

stopifnot(!anyDuplicated(
    candidate_support[c("pair_label", "scenario")]
))

# Explicitly distinguish untested combinations from negative results.
matrix_data <- candidate_support |>
    select(pair_label, scenario, support) |>
    mutate(support = as.character(support)) |>
    complete(
        pair_label,
        scenario,
        fill = list(support = "Not tested")
    ) |>
    mutate(
        pair_label = factor(pair_label, levels = rev(pair_order)),
        scenario = factor(scenario, levels = scenario_key$scenario),
        support = factor(
            support,
            levels = c(support_levels, "Not tested")
        )
    )

mediation_xenium_support_matrix <- ggplot(
    matrix_data,
    aes(x = scenario, y = pair_label)
) +
    geom_tile(aes(fill = support), color = "white", linewidth = 0.5) +
    geom_point(
        data = filter(matrix_data, support == "Not tested"),
        shape = 4, size = 1.5, color = "grey65"
    ) +
    scale_fill_manual(
        values = c(support_colors, "Not tested" = "white"),
        drop = FALSE
    ) +
    scale_x_discrete(
        labels = setNames(
            scenario_key$scenario_label,
            scenario_key$scenario
        )
    ) +
    labs(
        x = "Xenium scenario: mediator population → outcome population",
        y = NULL,
        fill = "Within-scenario support",
        title = "Xenium support for candidate mediation pairs",
        caption = paste(
            "Nominal P < .10. Scenarios reuse donors;",
            "multiple supported scenarios are not independent replications."
        )
    ) +
    theme_minimal(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom"
    ) +
    guides(fill = guide_legend(ncol = 1))

ggsave(
    here(plot_dir, "mediation_xenium_support_matrix.png"),
    mediation_xenium_support_matrix,
    width = 13, height = 8, dpi = 300
)

mediation_xenium_deg <- read_csv(here("processed-data",  "21_Xenium", "18_xenium_DEG_mediation_LR","Xenium_mediator_outcome_DEG_results.csv"))

mediation_xenium_tested_snakey <- mediation_xenium_deg |>
    mutate(outcome = ifelse(validate_outcome, paste0(outcome, "*"), outcome),
           mediator = ifelse(validate_mediator, paste0(mediator, "*"), mediator)) |>
    ggplot(aes(axis1=mediator, 
               axis2=outcome,
               fill=mediator
    )) +
    geom_alluvium() +
    geom_stratum(width=1/4) +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3) +
    # scale_fill_manual(values=c(
    #     "Astro.1"="#228B22",
    #     "Astro.2"="#32CD32",
    #     "Astro.3"="#808000"
    # )) +
    scale_x_discrete(limits=c("Mediator", "Outcome"), expand=c(0.1, 0.1)) +
    theme_void() +
    labs(fill="Astrocyte cluster", alpha="Validation status") +
    facet_wrap(~cluster_validate, nrow =1) +
    theme(legend.position = "None")

ggsave(mediation_xenium_tested_snakey, filename = here(plot_dir, "mediation_xenium_tested_snakey.png"), width = 9, height = 5)

mediation_xenium_hits <- read_csv(here("processed-data",  "22_Mediation", "03_Mediation_Xenium", "Xenium_mediation_hits.csv")) |> 
    mutate(pair = paste0(mediator, "|", outcome))


mediation_xenium_hits |> filter(mediator == "NPTXR", outcome == "GAD1") |> select(mediator, outcome, outcome_datatype, mediatorDE_P.Value) 
mediation_xenium_hits |> filter(mediator == "FZD8") |> select(mediator, outcome, med_cl_test, outcome_cl, mediatorDE_P.Value, mediated) |> arrange(outcome)

mediation_xenium_hits_snakey <- mediation_xenium_hits |>
    group_by(med_cl, mediator, outcome) |>
    summarize(n = n(),
              mediatorDE_P.Value_min = min(mediatorDE_P.Value, na.rm = TRUE),
              mediatorDE_P.Value_max = max(mediatorDE_P.Value, na.rm = TRUE),
              mediated = any(mediated)) |>
    mutate(
        outcome = ifelse(mediated, paste0(outcome, "*"), outcome),
           mediator = ifelse(mediatorDE_P.Value_min < 0.1, paste0(mediator, "*"), mediator)) |>
    ggplot(aes(axis1=mediator, 
               axis2=outcome,
               fill=mediator
    )) +
    geom_alluvium() +
    geom_stratum(width=1/2) +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), size=4) +
    scale_x_discrete(limits=c("Mediator", "Outcome"), expand=c(0.1, 0.1)) +
    theme_void() +
    labs(fill="Meditor", alpha="Validation status") +
    # facet_wrap(~med_cl, ncol =1)  +
    theme(legend.position = "None")

ggsave(mediation_xenium_hits_snakey, filename = here(plot_dir, "mediation_xenium_hits_snakey.png"), width = 3, height = 4)


#### ScDotPlots ####

message(Sys.time(), " - Load HDF5 sce")
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce_objects", "sce_ERC_subcluster"))
sce 

rownames(sce) <- rowData(sce)$gene_name

cell_type_colors <- metadata(sce)$cell_type_colors

all_mediator_genes <- unique(mediator_outcome$mediator)
all_outcome_genes <- unique(mediator_outcome$outcome)

sce$cell_type_mediation <- as.character(sce$cell_type_anno)

sce$cell_type_mediation[!sce$cell_type_broad %in% c("Astro", "Oligo")] <- as.character(sce$cell_type_broad[!sce$cell_type_broad %in% c("Astro", "Oligo")])
table(sce$cell_type_mediation)

pdf(here(plot_dir, sprintf("sn_cell_type_broad_dotplot_mediator_genes.pdf")), height = 5, width = 5)

sce|>
    scDotPlot(features = all_mediator_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = TRUE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = all_mediator_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = FALSE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()

pdf(here(plot_dir, sprintf("sn_cell_type_broad_dotplot_outcome_genes.pdf")), height =12, width = 5)

sce|>
    scDotPlot(features = all_outcome_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = TRUE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = all_outcome_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = FALSE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()


xenium_mediator_genes <- unique(mediation_xenium_deg$mediator)
xenium_outcome_genes <- unique(mediation_xenium_deg$outcome)


pdf(here(plot_dir, sprintf("sn_cell_type_broad_dotplot_mediator_genes_Xen.pdf")), height = 5, width = 5)

sce|>
    scDotPlot(features = xenium_mediator_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = TRUE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = xenium_mediator_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = FALSE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()

pdf(here(plot_dir, sprintf("sn_cell_type_med_dotplot_mediator_genes_Xen.pdf")), height = 5, width = 5)

sce|>
    scDotPlot(features = xenium_mediator_genes,
              group = "cell_type_mediation",
              groupAnno = "cell_type_mediation",
              scale = TRUE,
              # annoColors = list("cell_type_mediation" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = xenium_mediator_genes,
              group = "cell_type_mediation",
              groupAnno = "cell_type_mediation",
              scale = FALSE,
              # annoColors = list("cell_type_mediation" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()

pdf(here(plot_dir, sprintf("sn_cell_type_broad_dotplot_outcome_genes_Xen.pdf")), height = 9, width = 5)

sce|>
    scDotPlot(features = xenium_outcome_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = TRUE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = xenium_outcome_genes,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = FALSE,
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()


pdf(here(plot_dir, sprintf("sn_cell_type_med_dotplot_outcome_genes_Xen.pdf")), height = 9, width = 5)

sce|>
    scDotPlot(features = xenium_outcome_genes,
              group = "cell_type_mediation",
              groupAnno = "cell_type_mediation",
              scale = TRUE,
              # annoColors = list("cell_type_mediation" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = xenium_outcome_genes,
              group = "cell_type_mediation",
              groupAnno = "cell_type_mediation",
              scale = FALSE,
              # annoColors = list("cell_type_mediation" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()

### medation hits 

xenium_hit_gene_tb <- mediation_xenium_hits  |>
    distinct(pair, mediator, outcome) |>
    pivot_longer(!pair, names_to = "role", values_to = "gene") |>
    group_by(gene) |>
    summarise(role = ifelse(n()==2, "both", role)) |>
    add_row(gene = c("SOX5", "SOX6"), role = "outcome")

rowData(sce)$MedRole <- NULL
rowData(sce)$MedRole <- xenium_hit_gene_tb$role[match(rownames(sce), xenium_hit_gene_tb$gene)] 
table(rowData(sce)$MedRole)


pdf(here(plot_dir, sprintf("sn_cell_type_broad_dotplot_mediator_hit_Xen.pdf")), height = 4.5, width = 4.5)

sce|>
    scDotPlot(features = xenium_hit_gene_tb$gene,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = TRUE,
              featureAnno = "MedRole",
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

sce|>
    scDotPlot(features = xenium_hit_gene_tb$gene,
              group = "cell_type_broad",
              groupAnno = "cell_type_broad",
              scale = FALSE,
              featureAnno = "MedRole",
              annoColors = list("cell_type_broad" = cell_type_colors$broad),
              clusterRows = TRUE,
              groupLegends = FALSE
    )

dev.off()


