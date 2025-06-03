## Louise Huuki-Myers, May 2025
## Run  dreamlet on sn data

library("SpatialExperiment")
library("dreamlet")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "10_dreamlet_sn", "02_run_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "02_run_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- saveRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

#### Differential Expression ####
message(Sys.time(), ' - Differential Expression')

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round))

saveRDS(res.dl, file = here(data_dir, "dreamlet_sn.RDS"))

# names of estimated coefficients
coefNames(res.dl)

# results from full analysis
topTable(res.dl, coef = "APOE_carrierE4+")

pdf(here(plot_dir, "sn_dreamlet_VolcanoPlot.pdf"), height = 11, width = 8)
plotVolcano(res.dl, coef = "APOE_carrierE4+")
dev.off()

#### plot genes ####
# get data
df <- extractData(res.proc, "Inhib", genes = "NBPF12")

# expression boxplot
expression_plot <- ggplot(df, aes(APOE_carrier, NBPF12)) +
    geom_boxplot() +
    ylab(bquote(Expression ~ (log[2] ~ CPM))) +
    ggtitle("NBPF12") +
    theme_bw()

ggsave(expression_plot, filename = here(plot_dir, "sn_dreamlet_expression_boxplot.png"))

## forest plot
plotForest(res.dl, coef = "APOE_carrierE4+", gene = "NBPF12")

# slurmjobs::job_single('02_run_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 02_run_dreamlet_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

