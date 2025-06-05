## Louise Huuki-Myers, May 2025
## Run dreamlet Variance Parition on Visium data

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
data_dir <- here("processed-data", "11_dreamlet_Visium", "02_VarPart_dreamlet_Visium")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "11_dreamlet_Visium", "02_VarPart_dreamlet_Visium")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "11_dreamlet_Visium", "01_prep_dreamlet_Visium", "Visium_res_proc.rds"))

##### run variance partitioning ####
# source(here("code", "11_dreamlet_Visium", "dreamlet_models_Visium.R"))

dreamlet_models_Visium <- list(
    carrier = ~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide) + subsets_Mito_percent,
    apoe = ~ (1 | APOE) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | Visium_slide) + subsets_Mito_percent,
    e4e4 = ~ (1 | APOE_E4E4) + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | Visium_slide) + subsets_Mito_percent,
    # interaction syntax x + (x | g)
    carrier_i = ~ Anc_Afr + (Anc_Afr | APOE_carrier) + (1 | Sex) + Age + Rin + (1 | Visium_slide) + subsets_Mito_percent,
    apoe_i = ~ Anc_Afr + (Anc_Afr | APOE) + (1 | Sex) + Age + Rin + (1 | Visium_slide) + subsets_Mito_percent,
    e4e4_i = ~ Anc_Afr + (Anc_Afr | APOE_E4E4) + (1 | Sex) + Age + Rin + (1 | Visium_slide) + subsets_Mito_percent
)

stopifnot(opt$model %in% names(dreamlet_models_Visium))

message(Sys.time(), ' - Variance Partition, model = ', opt$model)
mod = dreamlet_models_Visium[[opt$model]]
print(mod)

vp.lst <- fitVarPart(res.proc, mod) 

saveRDS(vp.lst, file = here(data_dir, sprintf("dreamlet_fitVarPart_Visium-%s.RDS", opt$model)))

### VarPart plots ####
message(Sys.time(), " - Variance Partition plots")

AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    dplyr::filter(symbol %in% rownames(assay(res.proc, 1))) 

# Show variance fractions at the gene-level for each cell type
vp.lst[order(vp.lst[[1]], decreasing = TRUE)[1:100],]

coef_list <- list(carrier = "APOE_carrier",
                  carrier_i = "APOE_carrier.Anc_Afr",
                  apoe = "APOE",
                  apoe_i = "APOE.Anc_Afr",
                  e4e4 = "APOE_E4E4",
                  e4e4_i = "APOE_E4E4.Anc_Afr")

coef <- coef_list[[opt$model]]

top20_coef_varPart <- vp.lst$gene[order(vp.lst[[coef]], decreasing = TRUE)[1:20]]

pdf(here(plot_dir, sprintf("Visium_dreamlet_VarPart-%s.pdf", opt$model)), height = 11)
plotVarPart(vp.lst, label.angle = 45) + labs(title = "VarPart", subtitle = opt$model)
plotPercentBars(vp.lst[vp.lst$gene %in% AD_risk$symbol, ]) + labs(title = "Risk Genes", subtitle = opt$model)
plotPercentBars(vp.lst[vp.lst$gene %in% top20_coef_varPart, ]) + labs(title = "Top20 ", subtitle = opt$model)
dev.off()

# slurmjobs::job_single('02_VarPart_dreamlet_Visium', create_shell = TRUE, memory = '50G', command = "Rscript 02_VarPart_dreamlet_Visium.R")

slurmjobs::job_loop(
    loops = list(model = names(dreamlet_models_Visium)),
    name = "02_VarPart_dreamlet_Visium",
    create_shell = TRUE,
    create_script = FALSE
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

