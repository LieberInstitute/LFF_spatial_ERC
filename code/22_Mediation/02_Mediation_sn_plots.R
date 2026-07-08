## Louise Huuki-Myers, July 2026
## Compile an re-format mediation data - plot top pairs

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "22_Mediation", "02_Mediation_sn_plots")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "22_Mediation", "02_Mediation_sn_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### DEGs - Astro Mediation Genes ####

mediator_outcome <- read.delim(here("processed-data","22_Mediation","out-erc_astro","mediator_outcome_fdr_impact.tsv.gz")) |>
    as_tibble() |>
    mutate(pair = paste0(mediator, "|", outcome))

mediator_outcome |> count(fdr < 0.05, fdr_med > 0.05)
#      `fdr < 0.05` `fdr_med > 0.05`     n
#      <lgl>        <lgl>            <int>
#     1 TRUE         TRUE               619

# med_cl  mediator outcome     fdr logFC     t fdr_med logFC_med t_med pair          
# <chr>   <chr>    <chr>     <dbl> <dbl> <dbl>   <dbl>     <dbl> <dbl> <chr>         
# 1 Astro.1 NPTXR    ARHGAP31 0.0232 1.75   3.84  0.144      0.805  2.26 NPTXR|ARHGAP31
# 2 Astro.1 NPTXR    CALM3    0.0421 0.938  3.21  0.333      0.350  1.51 NPTXR|CALM3   
# 3 Astro.1 NPTXR    CCND2    0.0238 1.48   3.79  0.0938     0.806  2.66 NPTXR|CCND2   
# 4 Astro.1 NPTXR    CNR1     0.0284 1.76   3.59  0.0701     1.02   2.91 NPTXR|CNR1 

