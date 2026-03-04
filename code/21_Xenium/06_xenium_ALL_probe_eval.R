## Louise Huuki-Myers, March 2026
## Build custom list

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("readxl")
library("writexl")

data_dir <- here("processed-data", "21_Xenium", "06_xenium_ALL_probe_eval")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "06_xenium_ALL_probe_eval")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load probe data ####

ALL_probes_summary <- read_xlsx(here("processed-data", "21_Xenium", "02_xenium_compile_custom", "ERC_Xenium_ALL_probes_summary.xlsx")) |>
    filter(yesprobe) ## fix

select_custom_probes <- read_xlsx(here("processed-data", "21_Xenium", "02_xenium_compile_custom", "ERC_Xenium_select_custom_probes.xlsx"))

ALL_probes_summary |> count(in_base)
ALL_probes_summary |> filter(!in_base, !gene_name %in% select_custom_probes$gene_name)

ALL_probes_summary |> filter(grepl("SOX", gene_name))

#### DEGs - Astro Mediation Genes ####

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    mutate(mediator_probe = outcome %in% ALL_probes_summary$gene_name,
           outcome_probe = outcome %in% ALL_probes_summary$gene_name,
           pair = paste0(mediator, "|", outcome))

mediator_outcome_eval <- mediator_outcome |> filter(mediator_probe & outcome_probe)

mediator_outcome_eval |> unique()

write_csv(mediator_outcome_eval, file = here(data_dir, "ECR_mediator_outcome_xenium_eval.csv"))


mediator_outcome_eval |> count(mediator) # 9 mediators
mediator_outcome_eval |> count(outcome) # 31 outcome genes

ALL_probes_summary |>filter(gene_name %in% c(mediator_outcome_eval$mediator, mediator_outcome_eval$outcome)) |> count(goals)

