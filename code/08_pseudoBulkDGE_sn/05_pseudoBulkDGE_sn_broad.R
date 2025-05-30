## Louise Huuki-Myers, May 2025
## Run  pseudoBulkDGE on cell types from sn data

library("SingleCellExperiment")
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
sce_pb$APOE_E4E4 <- sce_pb$APOE == "E4/E4"

#### Design model ####

models <- list(carrier_n0 = ~APOE_carrier,
               carrier_n1 = ~APOE_carrier + ncells,
               carrier_n2 = ~APOE_carrier + pseudo_expr_chrM_ratio,
               carrier_n3 = ~APOE_carrier + Anc_Afr,
               carrier_n4 = ~APOE_carrier + ncells + pseudo_expr_chrM_ratio ,
               carrier = ~APOE_carrier + Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               carrier_i = ~APOE_carrier*Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               ## APOE
               APOE_n0 = ~APOE,
               APOE_n1 = ~APOE + ncells + pseudo_expr_chrM_ratio,
               APOE = ~APOE + Anc_Afr + Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               APOE_i = ~APOE*Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               APOE_ni0 = ~APOE*Anc_Afr,
               APOE_ni1 = ~APOE*Anc_Afr + ncells + pseudo_expr_chrM_ratio,
               ## E4E4
               E4E4 = ~APOE_E4E4 + Anc_Afr + Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round,
               E4E4_i = ~APOE_E4E4*Anc_Afr + Age + Sex + Rin + ncells + pseudo_expr_chrM_ratio + exp_round
)

# map(models, ~colnames(model.matrix(.x , colData(sce_pb))))
# colnames(model.matrix(models$E4E4 , colData(sce_pb)))

dge_design <- list(
    ## carrier
    carrier_n0 = list(mod = models$carrier_n0, coef = "APOE_carrierE4+"),
    carrier_n1 = list(mod = models$carrier_n1, coef = "APOE_carrierE4+"),
    carrier_n2 = list(mod = models$carrier_n2, coef = "APOE_carrierE4+"),
    carrier_n3 = list(mod = models$carrier_n3, coef = "APOE_carrierE4+"),
    carrier_n4 = list(mod = models$carrier_n4, coef = "APOE_carrierE4+"),
    carrier = list(mod = models$carrier, coef = "APOE_carrierE4+"),
    carrier_i = list(mod = models$carrier_i, coef = "APOE_carrierE4+:Anc_Afr"),
    ## APOE
    APOE_n = list(mod = models$APOE_n0, coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4")),
    APOE = list(mod = models$APOE, coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4")),
    APOE_i = list(mod = models$APOE_i, coef = c("APOEE2/E3:Anc_Afr","APOEE3/E4:Anc_Afr","APOEE4/E4:Anc_Afr")),
    APOE_ni0 = list(mod = models$APOE_ni0, coef = c("APOEE2/E3:Anc_Afr","APOEE3/E4:Anc_Afr","APOEE4/E4:Anc_Afr")),
    APOE_ni1 = list(mod = models$APOE_ni1, coef = c("APOEE2/E3:Anc_Afr","APOEE3/E4:Anc_Afr","APOEE4/E4:Anc_Afr")),
    ## E4E4 TODO update mod
    E4E4 = list(mod = models$E4E4, coef = "APOE_E4E4TRUE"),
    E4E4_i = list(mod = models$E4E4_i, coef = "APOE_E4E4TRUE:Anc_Afr"),
    ## Sex
    APOE_Sex = list(mod = models$APOE, coef = "SexM"),
    ## Anc
    APOE_Anc = list(mod = models$APOE, coef = "Anc_Afr")
)

# map(dge_design, "mod")

#### Run DGE ####
de.results <- map2(dge_design, names(dge_design), function(design, des_name){
    
    message(Sys.time(), " - Run pseudoBulkDGE on model: ", des_name , " coef: ", paste(design$coef, collapse = "+"))
    
    de.results <- pseudoBulkDGE(
        sce_pb,
        label = sce_pb$cell_type_broad,
        condition = sce_pb$APOE, # controls filtering, use consistently
        design = design$mod, # map model
        coef = design$coef, #map coef
        row.data = rowData(sce_pb),
        method = "voom"
    )
    
    return(de.results)
})

# de.results$test <- pseudoBulkDGE(
#     sce_pb,
#     label = sce_pb$cell_type_broad,
#     condition = sce_pb$APOE, # controls filtering, use consistently
#     design = models$APOE,
#     coef = c("APOEE2/E3","APOEE3/E4","APOEE4/E4"),
#     row.data = rowData(sce_pb),
#     method = "edgeR"
# )

ran <- map_int(de.results, length) > 0
if(!all(ran)) message("!! models with error: ", paste(names(ran)[!ran], collapse = ","))

(test_summary <- map(de.results[ran], ~summarizeTestsPerLabel(decideTestsPerLabel(.x, threshold=0.1))))

test_summary2 <- map2_dfr(test_summary, names(test_summary), function(ts, name){
    ts <- as.data.frame(ts)
    colnames(ts) <- paste0("deg",colnames(ts))
    
    if(!"deg-1" %in% colnames(ts)) ts$`deg-1` <- NA
    if(!"deg1" %in% colnames(ts)) ts$deg1 <- NA
    
    ts <- ts |> 
        rownames_to_column("cell_type_broad") |> 
        rowwise() |>
        mutate(mod = name,
               deg_total = sum(`deg-1`, deg1, na.rm = TRUE))
    
    return(ts)
})

test_summary2 |> group_by(mod) |> summarise(sum(deg_total))

write.csv(test_summary2, file = here(data_dir, "sn_pseudoBulkDGE_summary.csv"))

## save
saveRDS(de.results, file = here(data_dir, "sn_pseudoBulkDGE.rds"))

#### pairwise by APOE ####

apoe_combn <- as.data.frame(combn(levels(sce_pb$APOE), 2))

apoe_combn_sort <- map(apoe_combn, ~sort(.x))
names(apoe_combn_sort) <- map_chr(apoe_combn_sort, ~paste0(.x, collapse = "-"))

de.results.APOE_pairwise <- map(apoe_combn_sort, function(apoe_pair){
    
    # apoe_pair <- sort(apoe_pair)
    apoe_coef <- paste0("APOE", apoe_pair[[2]])
    
    spe_temp <- sce_pb[, sce_pb$APOE %in% apoe_pair]
    spe_temp$APOE <- droplevels(spe_temp$APOE)
    
    # check model matrix
    mm <- model.matrix(dge_design$APOE$mod , colData(spe_temp))
    stopifnot(apoe_coef == colnames(mm)[2])
    
    message(Sys.time(), sprintf(" - pairwise DE: %s(%i) vs. %s(%i) (n=%i), coef = %s", 
                                apoe_pair[[1]], sum(mm[,apoe_coef] == 0), 
                                apoe_pair[[2]], sum(mm[,apoe_coef] == 1),
                                ncol(spe_temp), apoe_coef))
    
    # table(spe_temp$APOE)
    # table(spe_temp$APOE, spe_temp$cell_type_broad)
    
    de.results <- pseudoBulkDGE(
        spe_temp,
        label = spe_temp$cell_type_broad,
        condition = spe_temp$APOE,
        design = dge_design$APOE$mod, # APOE
        coef = apoe_coef,
        row.data = rowData(spe_temp),
        method = "voom"
    )
    return(de.results)
})

ran <- map_int(de.results.APOE_pairwise, length) > 0
if(!all(ran)) message("!! models with error: ", paste(names(ran)[!ran], collapse = ", "))

# !! models with error: E2/E2-E2/E3, E2/E2-E4/E4, E2/E3-E4/E4

(pairwise_test_summary <- map(de.results.APOE_pairwise[ran], ~summarizeTestsPerLabel(decideTestsPerLabel(.x, threshold=0.1))))

pairwise_test_summary2 <- map2_dfr(pairwise_test_summary, names(pairwise_test_summary), function(ts, name){
    ts <- as.data.frame(ts)
    colnames(ts) <- paste0("deg",colnames(ts))
    
    if(!"deg-1" %in% colnames(ts)) ts$`deg-1` <- NA
    if(!"deg1" %in% colnames(ts)) ts$deg1 <- NA
    
    ts <- ts |> 
        rownames_to_column("cell_type_broad") |> 
        rowwise() |>
        mutate(mod = "APOE",
               pair = name,
               deg_total = sum(`deg-1`, deg1, na.rm = TRUE))
    
    return(ts)
})

pairwise_test_summary2 |> group_by(mod, pair) |> summarise(sum(deg_total))

write.csv(pairwise_test_summary2, file = here(data_dir, "sn_pseudoBulkDGE_pairwise_summary.csv"))

## save
saveRDS(de.results.APOE_pairwise, file = here(data_dir, "sn_pseudoBulkDGE_APOE_pairwise.rds"))

# slurmjobs::job_single('05_pseudoBulkDGE_sn_broad', create_shell = TRUE, memory = '50G', command = "Rscript 05_pseudoBulkDGE_sn_broad.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
