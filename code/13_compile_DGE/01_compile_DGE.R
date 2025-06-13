## Louise Huuki-Myers, June 2025
## Compile and plot all DGE data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "13_compile_DGE", "01_compile_DGE")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "13_compile_DGE", "01_compile_DGE")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

load(here("processed-data", "project_colors.Rdata"))
AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) 
main_mods <- c("carrier", "apoe", "e4e4", "carrier_i", "apoe_i", "e4e4_i")

#### pseudobulkDGE data ####

pseudobulkDGE_fn <- list(sn = here("processed-data", "08_pseudoBulkDGE_sn", "05_pseudoBulkDGE_sn_broad", "sn_pseudoBulkDGE.rds"),
                          Visium = here("processed-data", "09_pseudoBulkDGE_Visium", "03_pseudoBulkDGE_Visium", "Visium_pseudoBulkDGE.rds"))

map_lgl(pseudobulkDGE_fn, file.exists)

## load main mod data
pseudobulkDGE_data <- map(pseudobulkDGE_fn, ~readRDS(.x)[c("carrier", "APOE", "E4E4", "carrier_i", "APOE_i", "E4E4_i")])
## fix names
names(pseudobulkDGE_data$sn) <- main_mods
names(pseudobulkDGE_data$Visium) <- main_mods


pseudobulkDGE_data$sn$carrier$Astro

edit_pseudobulkDEG_data <- function(df, cluster){
    
    df2 <- df |>
        as.data.frame() |>
        mutate(cluster = cluster) |>
        select(gene_id, gene_name, logFC, t, P.Value. adj.P.Val)
    
    names_keep <- c()
    
}
    
    
#### dreamlet data ####
dreamlet_fn <- list(sn = here("processed-data", "10_dreamlet_sn", "05_compile_dreamlet_sn", "dreamlet_sn_TopTables.Rdata"),
                         Visium = here("processed-data", "11_dreamlet_Visium", "05_compile_dreamlet_Visium", "dreamlet_Visium_TopTables.Rdata"))

map_lgl(dreamlet_fn, file.exists)

dreamlet_data <- map(dreamlet_fn, ~get(load(.x)))

dreamlet_data$sn
