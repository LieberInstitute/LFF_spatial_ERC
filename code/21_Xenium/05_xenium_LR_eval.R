## Louise Huuki-Myers, Feb 2026
## Evaluate ligad receptor pairs in Xenium panel 

#### Set Up ####
# library(SingleCellExperiment)
library(tidyverse)
# library(liana)
library(here)
library(sessioninfo)
library(data.table)
library(readxl)

## Xenium probe set

ALL_probes_summary <- readxl::read_xlsx(here("processed-data", "21_Xenium", "02_xenium_compile_custom", "ERC_Xenium_ALL_probes_summary.xlsx"))

## liana data
#pre MFA non-spatially aware data
liana_data <- fread(here("processed-data", "17_cell_cell_comm", "liana", "ranked_results", "fine.csv.gz"))

liana_data_summary <- liana_data |>
    group_by(source, target, ligand_complex, receptor_complex) |>
    summarise(n_pass_magnitude_rank = sum(magnitude_rank < 0.05),
              n_pass_specificity_rank = sum(specificity_rank < 0.05, na.rm = TRUE),
              n_test = n()) 
# |>
#     mutate(source = factor(source, levels = cluster_levels),
#            target = factor(target, levels = cluster_levels))

liana_data_summary |> ungroup()|> filter(grepl("Astro", source) & target == "Oligo.3") |> arrange(-n_pass_magnitude_rank)

mediated_gene_shortlist <- c("MBP", # mylenation - ABCA8 outcome
                             "NAP1L3", # only outcome pf ST18
                             "ENC1", # only outcome of PLPPR4
                             "SOX6", # key regulator of oligodendrocyte maturation - NPTXR outcome
                             "SORBS1", # IL6ST outcome
                             "GPM6A" #Oligo.3 marker FZD8 outcome
)


mediated_gene_shortlist <- c("ADRB2")

liana_data_summary |> 
    ungroup() |>
    dplyr::filter(ligand_complex %in% mediated_gene_shortlist | ligand_complex %in% mediated_gene_shortlist) |>
    arrange(-n_pass_magnitude_rank) |>
    count(ligand_complex, receptor_complex)


# 1,779 interactions already covered
liana_xenium <- liana_data_summary |> 
    filter(ligand_complex %in% ALL_probes_summary$gene_name & 
               receptor_complex %in% ALL_probes_summary$gene_name)


liana_xenium |> ungroup()|> count(n_pass_magnitude_rank > 20) 
liana_xenium |> ungroup()|> count(n_pass_magnitude_rank) |> arrange(-n_pass_magnitude_rank)


liana_xenium |> ungroup()|> count(source, target) |> filter(target == "Oligo.3")

