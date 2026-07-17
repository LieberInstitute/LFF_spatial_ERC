## Louise Huuki-Myers, July 2026
## plot scdotplots, snakey, and upset plots for mediation genes 


#### Set up ####

library("tidyverse")
library("here")
library("ggalluvial")
library("scDotPlot")
library("SingleCellExperiment")
library("UpSetR")

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
              mediators = paste0(mediator, collapse = "|"))


outcome_summary |> filter(mediators == "CHRM3") |> dplyr::count(med_cl)
outcome_summary |> dplyr::count(mediators) ## 14 sets
    

#### Upset plots ####
# 
# mediated_groups <- mediator_outcome |>
#     mutate(
#         group = paste0(med_cl, "_", mediator)
#     ) |>
#     {\(df) split(df$outcome, df$group)}()
# 
# 
# pdf(here(plot_dir, "Mediation_gene_set_upset.pdf"))
# upset(fromList(mediated_groups), 
#       order.by = "freq", 
#       sets = names(mediated_groups), 
#       keep.order = TRUE
# )
# dev.off()

library(ComplexUpset)

mediator_outcome_wide <- mediator_outcome |> 
    mutate(mediator = paste0(med_cl, "_", mediator)) |>
    select(mediator, outcome) |>
    mutate(mediated = TRUE) |>
    pivot_wider(names_from = "mediator", values_from = "mediated") |>
    mutate(
        across(everything(), ~replace_na(.x, FALSE))
    )

mediator_outcome_wide |> filter(outcome == "SOX5")
mediator_outcome_wide |> filter(outcome %in% calcium_module) |> select(2:10)|> rowSums()

mediator_outcome_wide |> filter(Astro.2_ABCA8)

set_colors <- map(all_mediator_cl, ~upset_query(set=c(.x), fill=astro_colors[[jaffelab::ss(.x, "_")]]))

mediation_upsets <- ComplexUpset::upset(data = mediator_outcome_wide, 
                                        intersect = all_mediator_cl,
                                        name = "mediator",
                                        width_ratio=0.2,
                                        queries=set_colors) 

ggsave(mediation_upsets,filename =  here(plot_dir, "Mediation_gene_set_upset.pdf"), height = 4, width = 5)
ggsave(mediation_upsets,filename =  here(plot_dir, "Mediation_gene_set_upset.png"), height = 4, width = 6.5)


#### Snakey plots ####
mediator_snakey_all <- mediator_outcome |>
    ggplot(aes(x = mediator, next_node = outcome, next_x = outcome)) +
    geom_sankey()


# library("ggsankey")
# 
# mediation_xenium_deg_node <- mediation_xenium_deg |>
#     mutate(mediator_validate = ifelse(validate_mediator, "PASS", "FAIL"),
#            outcome_validate = ifelse(validate_outcome, "PASS", "FAIL")) |>
#     make_long(mediator, mediator_validate, outcome, outcome_validate)
# 
# 
# mediator_snakey_all <- mediation_xenium_deg_node |>
#     ggplot(aes(x = x, 
#                next_x = next_x, 
#                node = node, 
#                next_node = next_node,
#                fill = factor(node),
#                label = node)) +
#     geom_sankey() +
#     geom_sankey_label()

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

mediation_xenium_hits <- read_csv(here("processed-data",  "22_Mediation", "03_Mediation_Xenium", "Xenium_mediation_hits.csv"))


mediation_xenium_hits_snakey <- mediation_xenium_hits |>
    group_by(med_cl, mediator, outcome) |>
    summarize(n = n(),
              mediatorDE_P.Value_min = min(mediatorDE_P.Value),
              mediatorDE_P.Value_max = max(mediatorDE_P.Value),
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
    facet_wrap(~med_cl, nrow =1)  +
    theme(legend.position = "None")

ggsave(mediation_xenium_hits_snakey, filename = here(plot_dir, "mediation_xenium_hits_snakey.png"), width = 5, height = 4)


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

