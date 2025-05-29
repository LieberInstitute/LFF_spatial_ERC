## Louise Huuki-Myers, April 2025
## Run  pseudoBulkDGE

library("SpatialExperiment")
library("edgeR")
library("scran")
library("tidyverse")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("plots", "09_pseudoBulkDGE_Visium", "03_pseudoBulkDGE_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
spe_pb <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))

## temp
spe_pb$nspots <- spe_pb$ncells

#### Design model ####

models <- list(carrier = ~APOE_carrier + Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               carrier_i = ~APOE_carrier*Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               APOE = ~APOE + Anc_Afr + Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               APOE_i = ~APOE*Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide)

map(models, ~colnames(model.matrix(.x , colData(spe_pb))))


dge_design <- list(
    ## carrier
    carrier = list(mod = models$carrier, coef = "APOE_carrierE4+"),
    carrier_i = list(mod = models$carrier_i, coef = "APOE_carrierE4+:Anc_Afr"),
    ## APOE
    APOE = list(mod = models$APOE, coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4")),
    APOE_i = list(mod = models$APOE_i, coef = c("APOEE2/E3:Anc_Afr","APOEE3/E4:Anc_Afr","APOEE4/E4:Anc_Afr")),
    ## E4E4
    E4E4 = list(mod = models$APOE, coef = "APOEE4/E4"),
    E4E4_i = list(mod = models$APOE_i, coef = "APOEE4/E4:Anc_Afr")
)

#### Run DGE ####
de.results <- map2(dge_design, names(dge_design), function(design, des_name){
    
    message(Sys.time(), " - Run pseudoBulkDGE on model: ", des_name , " coef: ", design$coef)
    
    de.results <- pseudoBulkDGE(
        spe_pb,
        label = spe_pb$SpD,
        condition = spe_pb$APOE, # controls filtering, use consistently
        design = design$mod, # map model
        coef = design$coef, #map coef
        row.data = rowData(spe_pb),
        method = "edgeR"
    )
    
    return(de.results)
})

# de.results$test <- pseudoBulkDGE(
#     spe_pb,
#     label = spe_pb$SpD,
#     condition = spe_pb$APOE, # controls filtering, use consistently
#     design = models$APOE,
#     coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4"),
#     row.data = rowData(spe_pb),
#     method = "edgeR"
# )

map(de.results, ~summarizeTestsPerLabel(decideTestsPerLabel(.x, threshold=0.1)))

#### pariwise ####

apoe <- unique(spe_pb$APOE)
apoe_combn <- combn(apoe, 2)

map2(apoe_combn[1,], apoe_combn[2,], function(apoe1, apoe2){
    
    message(Sys.time(), sprintf("pairwise DE: %s vs. %s", apoe1, apoe2))
    spe_pb <- spe_pb[, spe_pb$APOE %in% c(apoe1, apoe2)]
    
    message(ncol(spe_pb))
    return(ncol(spe_pb))
    
})

### save ####
saveRDS(de.results, file = here(data_dir, "Visium_pseudoBulkDGE.rds"))

# slurmjobs::job_single('04_pseudoBulkDGE_Visium', create_shell = TRUE, memory = '50G', command = "Rscript 04_pseudoBulkDGE_Visium.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

