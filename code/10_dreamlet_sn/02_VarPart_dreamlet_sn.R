## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "02_VarPart_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "02_VarPart_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- saveRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

AD_risk <- read.csv(here("processed-data", "00_project_prep", "07_OpenTargets_AD_data", "clin_var_genes.csv")) |>
    dplyr::filter(symbol %in% rowData(sce)$gene_name) 

##### run variance partitioning ####
message(Sys.time(), ' - variance partition')

source(here("code", "10_dreamlet_sn", "dreamlet_models_sn.R"))
names(dreamlet_models_sn)

vp.lst <- fitVarPart(res.proc,~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + (1 | seq_round)) #subsets_Mito_percent + ncells

saveRDS(vp.lst, file = here(data_dir, "dreamlet_fitVarPart_sn.RDS"))
# vp.lst <- readRDS(here(data_dir, "dreamlet_fitVarPart_sn.RDS"))

# Show variance fractions at the gene-level for each cell type
# genes <- vp.lst$gene[2:4]

vp.lst[order(vp.lst$APOE_carrier, decreasing = TRUE)[1:100],]

top20_APOE_varPart <- vp.lst$gene[order(vp.lst$APOE_carrier, decreasing = TRUE)[1:20]]

pdf(here(plot_dir, "sn_dreamlet_VarPart.pdf"), height = 11)
plotVarPart(vp.lst, label.angle = 45)
plotPercentBars(vp.lst[vp.lst$gene %in% AD_risk$symbol, ]) + labs(title = "Risk Genes")
plotPercentBars(vp.lst[vp.lst$gene %in% top20_APOE_varPart, ]) + labs(title = "Top20 APOE")
dev.off()

# slurmjobs::job_single('02_VarPart_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 02_VarPart_dreamlet_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

