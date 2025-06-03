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
data_dir <- here("processed-data", "10_dreamlet_sn", "02_VarPart_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "02_VarPart_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

##### run variance partitioning ####
source(here("code", "10_dreamlet_sn", "dreamlet_models_sn.R"))
stopifnot(opt$model %in% names(dreamlet_models_sn))

message(Sys.time(), ' - Variance Partition, model = ', opt$model)
mod = dreamlet_models_sn[[opt$model]]
print(mod)

vp.lst <- fitVarPart(res.proc, mod) 

saveRDS(vp.lst, file = here(data_dir, sprintf("dreamlet_fitVarPart_sn-%s.RDS", opt$model)))

### VarPart plots ####
message(Sys.time(), " - Variance Partition plots")

AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    dplyr::filter(symbol %in% rownames(assay(res.proc, 1))) 

# Show variance fractions at the gene-level for each cell type
vp.lst[order(vp.lst[[1]], decreasing = TRUE)[1:100],]

top20_APOE_varPart <- vp.lst$gene[order(vp.lst[[1]], decreasing = TRUE)[1:20]]

pdf(here(plot_dir, sprintf("sn_dreamlet_VarPart-%s.pdf")), height = 11)
plotVarPart(vp.lst, label.angle = 45) + labs(title = "VarPart", subtitle = opt$model)
plotPercentBars(vp.lst[vp.lst$gene %in% AD_risk$symbol, ]) + labs(title = "Risk Genes", subtitle = opt$model)
plotPercentBars(vp.lst[vp.lst$gene %in% top20_APOE_varPart, ]) + labs(title = "Top20 ", subtitle = opt$model)
dev.off()

# slurmjobs::job_single('02_VarPart_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 02_VarPart_dreamlet_sn.R")

# slurmjobs::job_loop(
#     loops = list(model = names(dreamlet_models_sn)),
#     name = "02_VarPart_dreamlet_sn",
#     create_shell = TRUE,
#     create_script = FALSE
# )


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

