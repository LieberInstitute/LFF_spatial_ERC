## Louise Huuki-Myers, June 2026
## Explore Xenium DEG results for LR and mediation genes

#### Set Up ####
library("SpatialExperiment")
library("qs2")
library("here")
library("sessioninfo")
library("tidyverse")
library("spatialLIBD")
# library("DeconvoBuddies")
# library("scDotPlot")

data_dir <- here("processed-data", "21_Xenium", "18_xenium_DEG_mediation_LR")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "18_xenium_DEG_mediation_LR")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)



#### DEGs - Astro Mediation Genes ####

xenium_DEG <- readRDS(here("processed-data", "13_compile_DGE", "01_compile_DGE", "Xenium", "DGE_results_carrier_Xenium.Rds"))

# list.files(here("processed-data","22_Mediation","out-erc_astro"))

mediation_summary <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediation_vs_baseline_summary.tsv.gz"))

mediator_tab <- mediation_summary |> 
    filter(mediated_n > 0) |> 
    select(mediation_run, mediated_n) |> 
    separate(mediation_run, into = c("cell_type", "ensemblID", "gene_name"), sep = "\\|") |>
    mutate(in_base = gene_name %in% xenium_DEG$gene_name)

head(mediator_tab)

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(mediator_in_DEG = mediator %in% xenium_DEG$gene_name,
           outcome_in_DEG = outcome %in% xenium_DEG$gene_name,
           pair = paste0(mediator, "|", outcome))


xenium_DEG_med <- xenium_DEG |>
    inner_join(mediator_outcome |> select(cluster = med_cl, gene_name = mediator))


xenium_DEG_out <- xenium_DEG |> 
    filter(cluster == "Oligo.3") |>
    inner_join(mediator_outcome |> select(gene_name = outcome, pair, pair, fdr, logFC, t)) 

xenium_DEG_out |> filter(gene_name == "NPTXR") |>
    select(-c(6:12))

head(mediator_outcome)

