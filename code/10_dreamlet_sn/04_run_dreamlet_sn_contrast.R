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

ddf_suffix <- ifelse(opt$ddf == "Kenward-Roger", "_kr", "")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "03_run_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

#### Differential Expression ####

# The variable to be tested must be a fixed effect
dreamlet_contrast_models_sn <- list(
    contrast = ~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent
)
contrasts <-  c(
    `E2E2_E4E4` = "APOE_synE2.E2 - APOE_synE4.E4",
    `E3E4_E4E4` = "APOE_synE3.E4 - APOE_synE4.E4",
    `E2E3_E4E4` = "APOE_synE2.E3 - APOE_synE4.E4",
    `E2E2_E3E4` = "APOE_synE2.E2 - APOE_synE3.E4",
    `E2E2_E2E3` = "APOE_synE2.E2 - APOE_synE2.E3",
    `E2E3_E3E4` = "APOE_synE2.E3 - APOE_synE3.E4")



metadata(res.proc)

form <- ~ APOE + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium)

model.matrix(~0+APOE_syn, colData(res.proc))

# L <- makeContrastsDream(~0+APOE, colData(res.proc),
#                         contrasts = c(
#                             `E2/E2-E4/E4` = "APOEE2/E2 - APOEE4/E4",
#                             `E3/E4-E4/E4` = "APOEE3/E4 - APOEE4/E4",
#                             `E2/E3-E4/E4` = "APOEE2/E3 - APOEE4/E4",
#                             `E2/E2-E3/E4` = "APOEE2/E2 - APOEE3/E4",
#                             `E2/E2-E2/E3` = "APOEE2/E2 - APOEE2/E3",
#                             `E2/E3-E3/E4` = "APOEE2/E3 - APOEE3/E4")
# )

L <- makeContrastsDream(~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide), colData(res.proc),
                        contrasts = c(
                            `E2E2_E4E4` = "APOE_synE2.E2 - APOE_synE4.E4",
                            `E3E4_E4E4` = "APOE_synE3.E4 - APOE_synE4.E4",
                            `E2E3_E4E4` = "APOE_synE2.E3 - APOE_synE4.E4",
                            `E2E2_E3E4` = "APOE_synE2.E2 - APOE_synE3.E4",
                            `E2E2_E2E3` = "APOE_synE2.E2 - APOE_synE2.E3",
                            `E2E3_E3E4` = "APOE_synE2.E3 - APOE_synE3.E4")
)

# Visualize contrast matrix
contrast_plot <- plotContrasts(L) + labs(title = "DREAM contrast",
                                         subtitle = "~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide)")

ggsave(contrast_plot, filename = here(plot_dir, "Visium_dream_contrast.png"))



# source(here("code", "10_dreamlet_sn", "dreamlet_contrast_models_sn.R"))
stopifnot(opt$model %in% names(dreamlet_contrast_models_sn))

message(Sys.time(), ' - Differential Expression, model = ', opt$model, ", ddf = ", opt$ddf)
mod = dreamlet_contrast_models_sn[[opt$model]]
print(mod)

# Differential expression analysis within each assay,
# evaluated on the voom normalized data

param <- SnowParam(4, "SOCK", progressbar = TRUE)
res.dl <- dreamlet(res.proc, mod, ddf = opt$ddf, contrast)

message(Sys.time(), " - DONE Differential Expression...Save Data")
saveRDS(res.dl, file = here(data_dir, sprintf("dreamlet_sn-%s%s.RDS", opt$model, ddf_suffix)))

# names of estimated coefficients
message("coef")
coefNames(res.dl)

# topTable(res.dl, coef = "APOE_carrierE4+")

# slurmjobs::job_single('03_run_dreamlet_sn_kr', create_shell = TRUE, memory = '10G', command = "Rscript 03_run_dreamlet_sn.R --model carrier --ddf 'Kenward-Roger'")

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

