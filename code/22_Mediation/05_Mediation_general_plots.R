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

#### Load Data ####
mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(pair = paste0(mediator, "|", outcome))


mediator_outcome |> count(mediator)

#### Upset plots ####

mediated_groups <- mediator_outcome |>
    mutate(
        group = paste0(med_cl, "_", mediator)
    ) |>
    {\(df) split(df$outcome, df$group)}()


pdf(here(plot_dir, "Mediation_gene_set_upset.pdf"))
upset(fromList(mediated_groups), 
      order.by = "freq", 
      sets = names(mediated_groups), 
      keep.order = TRUE
)
dev.off()


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

