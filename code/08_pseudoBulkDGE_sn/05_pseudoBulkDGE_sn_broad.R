## Louise Huuki-Myers, May 2025
## Run  pseudoBulkDGE on cell types from sn data

library("SingleCellExperiment")
# library("edgeR")
library("scran")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "08_pseudoBulkDGE_sn", "05_pseudoBulkDGE_sn_broad")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
sce_pb <- readRDS(here("processed-data", "08_pseudoBulkDGE_sn", "01_pseudobulk_data_sn","sce_pseudo_DGE-cell_type_broad.RDS"))
dim(sce_pb)

table(sce_pb$cell_type_broad)

#### Design model + run DGE ####

models <- list(carrier = ~APOE_carrier + Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               carrier_i = ~APOE_carrier*Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               APOE = ~APOE + Anc_Afr + Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               APOE_i = ~APOE*Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round)

map(models, ~colnames(model.matrix(.x , colData(sce_pb))))


dge_design <- list(
    ## carrier
    carrier = c(mod = models$carrier, coef = "APOE_carrierE4+"),
    carrier_i = c(mod = models$carrier_i, coef = "APOE_carrierE4+:Anc_Afr"),
    ## APOE
    APOE = c(mod = models$APOE, coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4")),
    APOE_i = c(mod = models$APOE_i, coef = c("APOEE2/E3:Anc_Afr","APOEE3/E4:Anc_Afr","APOEE4/E4:Anc_Afr")),
    ## E4E4
    E4E4 = c(mod = models$APOE, coef = "APOEE4/E4"),
    E4E4_i = c(mod = models$APOE_i, coef = "APOEE4/E4:Anc_Afr")
)

de.results <- map2(dge_design, names(dge_design), function(design, des_name){
    
    message(Sys.time(), " - Run pseudoBulkDGE on model:", mod_name , " coef:", coef[[mod_name]])
    
    de.results <- pseudoBulkDGE(
        sce_pb,
        label = sce_pb$cell_type_broad,
        condition = sce_pb$carrier, # controls filtering, use consistently
        design = design$mod, # map model
        coef = dge_design$coef, #map coef
        row.data = rowData(sce_pb),
        method = "edgeR"
    )

    return(de.results)
})

# de.results$test <- pseudoBulkDGE(
#     sce_pb,
#     label = sce_pb$cell_type_broad,
#     condition = sce_pb$carrier, # controls filtering, use consistently
#     design = models$APOE,
#     coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4"),
#     row.data = rowData(sce_pb),
#     method = "edgeR"
# )

map(de.results, ~summarizeTestsPerLabel(decideTestsPerLabel(.x, threshold=0.1)))

saveRDS(de.results, file = here(data_dir, "sn_pseudoBulkDGE_broad.rds"))

# slurmjobs::job_single('05_pseudoBulkDGE_sn_broad', create_shell = TRUE, memory = '50G', command = "Rscript 05_pseudoBulkDGE_sn_broad.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
