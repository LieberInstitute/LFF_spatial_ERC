## Louise Huuki-Myers, July 2026
## Compile an re-format mediation data - plot top pairs

## set up
library("tidyverse")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "22_Mediation", "02_Mediation_sn_plots")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "22_Mediation", "02_Mediation_sn_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## DEGs - astro mediation genes

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    as_tibble() |>
    mutate(pair = paste0(mediator, "|", outcome))

mediator_outcome |> count(fdr_base < 0.05, fdr_med > 0.05, fdr_med_vec < 0.05)
## `*_base` are baseline APOE stats, `*_med` are APOE stats after adding
## `med_vec`, and `*_med_vec` are the mediator coefficient stats.
