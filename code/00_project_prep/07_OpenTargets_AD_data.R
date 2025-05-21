## Louise Huuki-Myers, May 2025
## Reformat and examine OpenTargets data

library("tidyverse")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "00_project_prep", "07_OpenTargets_AD_data")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# plot_dir <- here("plots", "00_project_prep", "07_OpenTargets_AD_data")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
open_target_tb <- read_tsv(here("external-data", "OpenTargets", "OT-MONDO_0004975-associated-targets-5_21_2025-v25_03.tsv")) 

(gene_sets <- colnames(open_target_tb)[-1])
# [1] "globalScore"                        "gwasCredibleSets"                   "geneBurden"                        
# [4] "eva"                                "genomicsEngland"                    "gene2Phenotype"                    
# [7] "uniprotLiterature"                  "uniprotVariants"                    "orphanet"                          
# [10] "clingen"                            "cancerGeneCensus"                   "intogen"                           
# [13] "evaSomatic"                         "cancerBiomarkers"                   "chembl"                            
# [16] "crisprScreen"                       "crispr"                             "slapenrich"                        
# [19] "progeny"                            "reactome"                           "sysbio"                            
# [22] "europepmc"                          "expressionAtlas"                    "impc"                              
# [25] "maxClinicalTrialPhase"              "isInMembrane"                       "isSecreted"                        
# [28] "hasLigand"                          "hasSmallMoleculeBinder"             "hasPocket"                         
# [31] "mouseOrthologMaxIdentityPercentage" "hasHighQualityChemicalProbes"       "geneticConstraint"                 
# [34] "mouseKoScore"                       "geneEssentiality"                   "hasSafetyEvent"                    
# [37] "isCancerDriverGene"                 "paralogMaxIdentityPercentage"       "tissueSpecificity"                 
# [40] "tissueDistribution"

select_gene_sets <- c("geneBurden", "eva", "genomicsEngland", "gene2Phenotype","uniprotLiterature","uniprotVariants","orphanet","clingen")

open_target_tb <- open_target_tb |> mutate_at(gene_sets[-1], ~as.double(na_if(., "No data")))

summary(open_target_tb)

## eva = ClinVar (?)

open_target_tb |> filter(symbol %in% c("PSEN1", "APP", "PSEN2", "HFE"))

clin_var_genes <- open_target_tb |> 
    filter(!is.na(eva) & eva > 0) |> 
    select(symbol, eva) 

clin_var_genes |> arrange(-eva) |> print(n = 25)

write_csv(clin_var_genes, file = here(data_dir, "clin_var_genes.csv"))

## make long and examine scores
open_target_tb_long <- open_target_tb |> 
    pivot_longer(!symbol, names_to = "gene_score", values_to = "score") |> 
    filter(!is.na(score))

open_target_tb_long |> group_by(gene_score)  |> count() |> arrange(-n)

open_target_tb_long |> 
    group_by(gene_score) |> 
    filter(gene_score %in% select_gene_sets)  |> 
    count() |> 
    arrange(-n)

## output select sets
open_target_tb_long |> 
    filter(gene_score %in% select_gene_sets) |>
    write_csv(file = here(data_dir, "OpenTargets_select_scores.csv"))

## bar plot
open_target_tb_long |> 
    ggplot(aes(x = fct_infreq(gene_score))) +
    geom_bar()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

open_target_tb_long |> 
    filter(gene_score %in% select_gene_sets)  |> 
    ggplot(aes(x = fct_infreq(gene_score))) +
    geom_bar()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

## score density plot
open_target_tb_long |> 
    ggplot(aes(x = score)) + 
    geom_density() + 
    facet_wrap(~gene_score)

open_target_tb_long |> 
    filter(gene_score %in% select_gene_sets)  |> 
    ggplot(aes(x = score)) + 
    geom_histogram(bins = 10) + 
    facet_wrap(~gene_score)



