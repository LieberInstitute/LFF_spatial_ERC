## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("model", "m", "1", "character", "Model name"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)
print(opt)

# test 
# opt$model <- "carrier"

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "03_run_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

#### Differential Expression ####

# The variable to be tested must be a fixed effect
dreamlet_models_sn <- list(
    carrier = ~ APOE_carrier + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    apoe = ~ APOE + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    e4e4 = ~ APOE_E4E4 + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round) ,
    contrast = ~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    # interaction syntax x + (x | g)
    carrier_i = ~ APOE_carrier*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    apoe_i = ~ APOE*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    e4e4_i = ~ APOE_E4E4*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round)
)

# source(here("code", "10_dreamlet_sn", "dreamlet_models_sn.R"))
stopifnot(opt$model %in% names(dreamlet_models_sn))

message(Sys.time(), ' - Differential Expression, model = ', opt$model)
mod = dreamlet_models_sn[[opt$model]]
print(mod)

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, mod)

message(Sys.time(), ' - DONE Differential Expression')

# names of estimated coefficients
message("coef")
coefNames(res.dl)

topTable(res.dl, coef = "APOE_carrierE4+")

message(Sys.time(), " - Save Data")
saveRDS(res.dl, file = here(data_dir, sprintf("dreamlet_sn-%s.RDS", opt$model)))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

