## Louise Huuki-Myers, April 2025
## Run  pseudoBulkDGE on Visium data

library("SpatialExperiment")
library("edgeR")
library("scran")
library("tidyverse")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "09_pseudoBulkDGE_Visium", "03_pseudoBulkDGE_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
spe_pb <- readRDS(here("processed-data", "09_pseudoBulkDGE_Visium", "01_pseudobulk_data_Visium", "spe_pseudo_DGE.RDS"))

## temp
# spe_pb$APOE_E4E4 <- spe_pb$APOE == "E4/E4"

#### Design model ####

models <- list(carrier_n0 = ~APOE_carrier,
               carrier_n1 = ~APOE_carrier + nspots,
               carrier_n2 = ~APOE_carrier + pseudo_expr_chrM_ratio,
               carrier_n3 = ~APOE_carrier + Anc_Afr,
               carrier_n4 = ~APOE_carrier + nspots + pseudo_expr_chrM_ratio ,
               carrier = ~APOE_carrier + Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               carrier_i = ~APOE_carrier*Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               ## APOE
               APOE_n0 = ~APOE,
               APOE_n1 = ~APOE + nspots + pseudo_expr_chrM_ratio,
               APOE = ~APOE + Anc_Afr + Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               APOE_i = ~APOE*Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               APOE_ni0 = ~APOE*Anc_Afr,
               APOE_ni1 = ~APOE*Anc_Afr + nspots + pseudo_expr_chrM_ratio,
               ## E4E4
               E4E4 = ~APOE_E4E4 + Anc_Afr + Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide,
               E4E4_i = ~APOE_E4E4*Anc_Afr + Age + Sex + Rin + nspots + pseudo_expr_chrM_ratio + Visium_slide
               )

# map(models, ~colnames(model.matrix(.x , colData(spe_pb))))
colnames(model.matrix(models$E4E4 , colData(spe_pb)))

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
    E4E4 = list(mod = models$APOE, coef = "APOE_E4E4TRUE"),
    E4E4_i = list(mod = models$APOE_i, coef = "APOE_E4E4TRUE:Anc_Afr"),
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
        spe_pb,
        label = spe_pb$SpD,
        condition = spe_pb$APOE, # controls filtering, use consistently
        design = design$mod, # map model
        coef = design$coef, #map coef
        row.data = rowData(spe_pb),
        method = "voom"
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

ran <- map_int(de.results, length) > 0
if(!all(ran)) message("!! models with error: ", paste(names(ran)[!ran], collapse = ","))
    
(test_summary <- map(de.results[ran], ~summarizeTestsPerLabel(decideTestsPerLabel(.x, threshold=0.1))))

test_summary2 <- map2_dfr(test_summary, names(test_summary), function(ts, name){
    ts <- as.data.frame(ts)
    colnames(ts) <- paste0("deg",colnames(ts))
    
    if(!"deg-1" %in% colnames(ts)) ts$`deg-1` <- NA
    if(!"deg1" %in% colnames(ts)) ts$deg1 <- NA
    
    ts <- ts |> 
        rownames_to_column("SpD") |> 
        rowwise() |>
        mutate(mod = name,
               deg_total = sum(`deg-1`, deg1, na.rm = TRUE))
    
    return(ts)
    })

test_summary2 |> group_by(mod) |> summarise(sum(deg_total))

write.csv(test_summary2, file = here(data_dir, "Visium_pseudoBulkDGE_summary.csv"))

## save
saveRDS(de.results, file = here(data_dir, "Visium_pseudoBulkDGE.rds"))


de.results$carrier$`L2.3~Sp09D01`

# de.results$APOE$`L2.3~Sp09D01`[de.results$APOE$`L2.3~Sp09D01`$FDR < 0.1,] 


#### pairwise ####

apoe_combn <- as.data.frame(combn(unique(spe_pb$APOE), 2))

apoe_combn_sort <- map(apoe_combn, ~sort(.x))
names(apoe_combn_sort) <- map_chr(apoe_combn_sort, ~paste0(.x, collapse = "-"))


de.results.APOE_pairwise <- map(apoe_combn_sort, function(apoe_pair){
    
    # apoe_pair <- sort(apoe_pair)
    apoe_coef <- paste0("APOE", apoe_pair[[2]])
    
    spe_temp <- spe_pb[, spe_pb$APOE %in% apoe_pair]
    
    # check model matrix
    mm <- model.matrix(dge_design$APOE$mod , colData(spe_temp))
    stopifnot(apoe_coef == colnames(mm)[2])
    
    message(Sys.time(), sprintf(" - pairwise DE: %s(%i) vs. %s(%i) (n=%i), coef = %s", 
                                apoe_pair[[1]], sum(mm[,apoe_coef] == 0), 
                                apoe_pair[[2]], sum(mm[,apoe_coef] == 1),
                                ncol(spe_temp), apoe_coef))
    
    # table(spe_temp$APOE)
    # table(spe_temp$APOE, spe_temp$SpD)
    
    de.results <- pseudoBulkDGE(
        spe_temp,
        label = spe_temp$SpD,
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

# !! models with error: 
# E2/E2-E4/E4, 
# E2/E3-E4/E4, 
# E2/E2-E2/E3

# 2025-05-29 14:54:37.408092 - pairwise DE: E2/E2 vs. E4/E4, coef = APOEE4/E4 ! "APOEE4/E4"
# 2025-05-29 14:54:37.492065 - pairwise DE: E3/E4 vs. E4/E4, coef = APOEE4/E4   "APOEE4/E4"
# 2025-05-29 14:54:37.554442 - pairwise DE: E2/E3 vs. E4/E4, coef = APOEE4/E4 ! "APOEE4/E4"
# 2025-05-29 14:54:37.622258 - pairwise DE: E2/E2 vs. E3/E4, coef = APOEE3/E4   "APOEE3/E4"
# 2025-05-29 14:54:37.68476 - pairwise DE: E2/E2 vs. E2/E3, coef = APOEE2/E3 !  "APOEE2/E3"
# 2025-05-29 14:54:37.741843 - pairwise DE: E2/E3 vs. E3/E4, coef = APOEE3/E4   "APOEE3/E4"

map_int(de.results.APOE_pairwise, length)

map(de.results.APOE_pairwise[ran], ~summarizeTestsPerLabel(decideTestsPerLabel(.x, threshold=0.1)))

## save
saveRDS(de.results.APOE_pairwise, file = here(data_dir, "Visium_pseudoBulkDGE_APOE_pairwise.rds"))


# slurmjobs::job_single('04_pseudoBulkDGE_Visium', create_shell = TRUE, memory = '50G', command = "Rscript 04_pseudoBulkDGE_Visium.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

