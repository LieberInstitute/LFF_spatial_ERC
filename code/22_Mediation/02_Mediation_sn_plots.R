## Louise Huuki-Myers, July 2026
## Compile and re-format mediation data - plot top pairs

## set up
library("tidyverse")
library("here")
library("sessioninfo")
library("edgeR")
library("limma")
library("SingleCellExperiment")
library("patchwork")
library("DeconvoBuddies")

source(here("code", "utils", "plot_DEG_express.R"))
source(here("code", "22_Mediation", "plot_DEG_mediated_express.R"))
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

data_dir <- here("processed-data", "22_Mediation", "02_Mediation_sn_plots")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "22_Mediation", "02_Mediation_sn_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### Load mediation data ####
mediator_outcome <- read.delim(here("processed-data", "22_Mediation", "out-erc_astro", "mediator_outcome_fdr_impact.tsv.gz")) |>
    as_tibble() |>
    mutate(pair = paste0(mediator, "|", outcome))

xenium_mediation <- read_csv(here("processed-data", "22_Mediation", "03_Mediation_Xenium", "mediation_summary-all_scenarios_Pval0.10.csv")) |>
    filter(Xen_valid, !Xen_carrier_sig)  |>
    mutate(pair = paste0(mediator, "|", outcome)) 

mediation_pairs <- list(
    c1 = c("ABCA8|MBP", "IL6ST|THBS1", "IL6ST|STAT4", "FZD8|IGFBP5"),
    c2 = c("CHRM3|ERMN", "NETO1|MGST1", "NETO1|NRG3", "SV2B|ERBB3", "SV2B|MAG"),
    xen = xenium_mediation |>
        pull(pair) |>
        unique()
)


mediator_outcome_select <- mediator_outcome |>
    mutate(
        outcome_cl  = outcome_cl,
        med_cl_test = med_cl,
        fdr_SN = fdr,        
        t_SN = t,
        fdr_SN_carrier = fdr_med, 
        t_SN_carrier = t_med,
        base_sig    = fdr_SN < 0.05,
        carrier_sig = fdr_SN_carrier < 0.05,
        SN_mediated = !carrier_sig
    ) |>
    filter(pair %in% unlist(mediation_pairs))



#### load data: sce_pb + DE_data ####

sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn", "sce_pseudo_DGE-cell_type_anno.RDS"))
rownames(sce_pb) <- rowData(sce_pb)$gene_name

## load DE data 
DE_data <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "sn_fine", "DGE_results_carrier_sn_fine.Rds"))

## DE model params -- matches 04_DEG_boxplots.R and 01_mediation_screening.R's baseline_formula_str ####
outcome_cl <- "Oligo.3"
batch <- "exp_round"
sn_mod <- ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio
sn_cleanY_P <- 4  ## keeps all 4 APOE_syn genotype-level columns


## test plot_DEG_express w/ defined variables
MBP_DEG_plot <- plot_DEG_express(sce = sce_pb,
                 stats = DE_data,
                 clus = outcome_cl,
                 gene = "MBP",
                 pval_col = "vlmf_adj.P.Val",
                 fc_col = "vlmf_logFC",
                 gene_col = "gene_name",
                 cluster_col = "registration_variable",
                 category_col = "APOE_carrier",
                 mod = sn_mod,
                 color_pal = APOE_carrier_colors,
                 plot_points = TRUE,
                 ncol = 2,
                 cleanY_P = 4)

ggsave(MBP_DEG_plot, filename = here(plot_dir, "test_Oligo.3_MBP_DEG_plot.png"))



# write_csv(mediation_pairs_stats, here(data_dir, "mediation_pairs_stats_sn.csv"))


#### Test mediation plots ####

med_plot_test <- plot_DEG_mediated_express(
    sce = sce_pb,
    sce_mediator = sce_pb,
    stats = mediation_pairs_stats,
    clus = "Oligo.3",
    med_clus = "Astro.2",
    mediator_gene = "ABCA8",
    gene = "MBP",
    category_col = "APOE_carrier",
    mod = ~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio,
    cleanY_P = 4,
    color_pal = APOE_carrier_colors,
    plot_points = TRUE,
    plot_mediator_panel = TRUE,
    mediator_stats = DE_data,
    anno_stat_suffix = "SN",
    signif_stat = "fdr"
)

ggsave(med_plot_test, filename = here(plot_dir, "Mediation_plot_test_MBP_ABCA8.png"), width = 8, height =5)

#### Loop over each (med_cl, mediator) group and plot ####

mediator_outcome_select2 <- mediator_outcome_select |> select(med_cl, mediator, outcome)

pdf(here(plot_dir, "mediation_boxplots_sn_pairs.pdf"), width = 10, height = 4)

pwalk(mediator_outcome_select2, function(med_cl, mediator, outcome) {

    message(Sys.time(), sprintf(" - plotting mediator=%s (%s) -> %s", mediator, med_cl, outcome))

    p <- tryCatch(
        plot_DEG_mediated_express(
            sce = sce_pb,
            sce_mediator = sce_pb,
            stats = mediator_outcome_select,
            clus = outcome_cl,
            med_clus = med_cl,
            mediator_gene = mediator,
            gene = outcome,
            mod = sn_mod,
            cleanY_P = sn_cleanY_P,
            color_pal = APOE_carrier_colors,
            plot_points = TRUE,
            plot_mediator_panel = TRUE,
            mediator_stats = DE_data,
            anno_stat_suffix = "SN",
            signif_stat = "fdr"
        ),
        error = function(e) {
            warning(sprintf("Failed to plot mediator=%s med_cl=%s: %s", mediator, med_cl, conditionMessage(e)))
            NULL
        }
    )

    if (!is.null(p)) {
        print(p + patchwork::plot_annotation(title = sprintf("%s (%s) -> %s", mediator, med_cl, outcome_cl)))
        ggsave(p, filename = here(plot_dir, sprintf(
            "mediation_sn_%s_%s_Oligo3_%s.png", gsub("\\.", "", med_cl), mediator, outcome
        )), height = 4, width = 8)
    }
})

dev.off()

message(Sys.time(), sprintf(" - Done: %d groups plotted covering %d pairs", nrow(plot_groups), nrow(mediation_pairs_stats)))

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
