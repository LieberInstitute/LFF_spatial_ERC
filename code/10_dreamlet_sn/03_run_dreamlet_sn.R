## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("model", "m", "1", "character", "Model name",
      "ddf", "d", "1", "character", "degress of freedom method"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

# test 
# opt$model <- "carrier"
# opt$ddf <- "Kenward-Roger"

## use default "Satterthwaite" for df approximation
if(is.null(opt$ddf)) opt$ddf <- "Satterthwaite"

print(opt)

ddf_suffix <- ifelse(opt == "Kenward-Roger", "_kr", "")

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
    carrier = ~ APOE_carrier + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round)  + subsets_Mito_percent,
    apoe = ~ APOE + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    e4e4 = ~ APOE_E4E4 + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    #r2
    carrier_n0 = ~ APOE_carrier,
    carrier_n1 = ~ APOE_carrier + Anc_Afr,
    carrier_n2 = ~ APOE_carrier + Anc_Afr + (1 | Sex),
    carrier_n3 = ~ APOE_carrier + Anc_Afr + (1 | Sex) + Age ,
    carrier_n4 = ~ APOE_carrier + (1 | exp_round)  + subsets_Mito_percent,
    carrier_sf = ~ APOE_carrier + Anc_Afr + Sex + Age + Rin + (1 | exp_round)  + subsets_Mito_percent,
    contrast = ~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    apoe_n0 = ~ APOE,
    e4e4_n0 = ~ APOE_E4E4,
    # interaction syntax x + (x | g)
    carrier_i = ~ APOE_carrier*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    apoe_i = ~ APOE*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    e4e4_i = ~ APOE_E4E4*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent
)

# source(here("code", "10_dreamlet_sn", "dreamlet_models_sn.R"))
stopifnot(opt$model %in% names(dreamlet_models_sn))

message(Sys.time(), ' - Differential Expression, model = ', opt$model, ", ddf = ", opt$ddf)
mod = dreamlet_models_sn[[opt$model]]
print(mod)

# Differential expression analysis within each assay,
# evaluated on the voom normalized data

param <- SnowParam(4, "SOCK", progressbar = TRUE)
res.dl <- dreamlet(res.proc, mod, ddf = opt$ddf, BPPARAM = param)

message(Sys.time(), " - DONE Differential Expression...Save Data")
saveRDS(res.dl, file = here(data_dir, sprintf("dreamlet_sn-%s%s.RDS", opt$model, ddf_suffix)))

# names of estimated coefficients
message("coef")
coefNames(res.dl)

# topTable(res.dl, coef = "APOE_carrierE4+")

slurmjobs::job_single('03_run_dreamlet_sn_kr', create_shell = TRUE, memory = '10G', command = "Rscript 03_run_dreamlet_sn.R --model carrier --ddf 'Kenward-Roger'")

# slurmjobs::job_loop(
#     loops = list(model = c("carrier_n0","carrier_n1","carrier_n2","carrier_n3","carrier_n4","carrier_sf","apoe_n0","e4e4_n0", "contrast")),
#     name = "03_run_dreamlet_sn_r2",
#     create_shell = TRUE,
#     create_script = FALSE
# )

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

