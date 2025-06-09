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
# opt$model <- "contrast"
# opt$ddf <- "Kenward-Roger"

## use default "Satterthwaite" for df approximation
if(is.null(opt$ddf)) opt$ddf <- "Satterthwaite"

print(opt)

ddf_suffix <- ifelse(opt$ddf == "Kenward-Roger", "_kr", "")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "04_run_dreamlet_contrast_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "04_run_dreamlet_contrast_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
pb <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_dreamlet_pb.rds"))

#### Differential Expression ####

# The variable to be tested must be a fixed effect
dreamlet_contrast_models_sn <- list(
    contrast = ~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent
)

contrasts <-  c(
    E2E2_E4E4 = "APOE_synE2.E2 - APOE_synE4.E4",
    E3E4_E4E4 = "APOE_synE3.E4 - APOE_synE4.E4",
    E2E3_E4E4 = "APOE_synE2.E3 - APOE_synE4.E4",
    E2E2_E3E4 = "APOE_synE2.E2 - APOE_synE3.E4",
    E2E2_E2E3 = "APOE_synE2.E2 - APOE_synE2.E3",
    E2E3_E3E4 = "APOE_synE2.E3 - APOE_synE3.E4",
    E4E4_anyE2 = "APOE_synE4.E4 - (0.5*APOE_synE2.E3 + 0.5*APOE_synE2.E2)",
    anyE4_anyE2 = "APOE_synE3.E4 - (0.5*APOE_synE2.E3 + 0.5*APOE_synE2.E2)",
    E2E2_anyE4 = "APOE_synE2.E2 - (0.5*APOE_synE4.E4 + 0.5*APOE_synE3.E4)",
    E2E3_anyE4 = "APOE_synE2.E3 - (0.5*APOE_synE4.E4 + 0.5*APOE_synE3.E4)")


stopifnot(opt$model %in% names(dreamlet_contrast_models_sn))

message('model = ', opt$model, ", ddf = ", opt$ddf)
mod = dreamlet_contrast_models_sn[[opt$model]]
print(mod)

message(Sys.time(), " - Apply Voom")
# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, mod, min.count = 5, min.cells = 10)

# show voom plot for each cell clusters
pdf(here(plot_dir, sprintf("sn_dreamlet_voom-%s.pdf", opt$model)))
plotVoom(res.proc)
dev.off()

model.matrix(~0+APOE_syn, colData(res.proc))

#### Visualize contrast matrix ####
L <- makeContrastsDream(~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round), 
                        colData(res.proc),
                        contrasts = contrasts
)

contrast_plot <- plotContrasts(L) + labs(title = "sn dreamlet contrast",
                                         subtitle = mod)

ggsave(contrast_plot, filename = here(plot_dir, "sn_dreamlet_contrast.png"))

message(Sys.time(), " - dreamlet")

# Differential expression analysis within each assay,
# evaluated on the voom normalized data

# param <- SnowParam(4, "SOCK", progressbar = TRUE)

res.dl <- dreamlet(res.proc, 
                   formula = ~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round), 
                   contrasts = contrasts)

## how to add mito rate to colData?

message(Sys.time(), " - DONE Differential Expression...Save Data")
saveRDS(res.dl, file = here(data_dir, sprintf("dreamlet_contrast_sn-%s%s.RDS", opt$model, ddf_suffix)))

# names of estimated coefficients
message("coef")
coefNames(res.dl)

# topTable(res.dl, coef = "APOE_carrierE4+")

# slurmjobs::job_single('04_run_dreamlet_contrast_sn', create_shell = TRUE, memory = '10G', command = "Rscript 04_run_dreamlet_contrast_sn.R --model contrast")

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

