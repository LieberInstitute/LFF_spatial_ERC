## Louise Huuki-Myers, June 2025
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
data_dir <- here("processed-data", "10_dreamlet_sn", "04_plot_dreamlet_sn")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "10_dreamlet_sn", "04_plot_dreamlet_sn")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### Load data ####
res.proc <- readRDS(here("processed-data", "10_dreamlet_sn", "01_prep_dreamlet_sn", "sn_res_proc.rds"))

res.dl <- readRDS(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn", sprintf("dreamlet_sn-%s.RDS", opt$model)))

res.dl.fn <- list.files(here("processed-data", "10_dreamlet_sn", "03_run_dreamlet_sn"), full.names = TRUE)
names(res.dl.fn) <- gsub(".RDS", "", gsub("dreamlet_sn-", "", basename(res.dl.fn)))

res.dl.list <- purrr::map(res.dl.fn, readRDS)
purrr::map(res.dl.list, coefNames)


details(res.dl)
coefNames(res.dl)

# results from full analysis
topTable(res.dl, coef = "Anc_Afr")

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

# slurmjobs::job_single('04_plot_dreamlet_sn', create_shell = TRUE, memory = '50G', command = "Rscript 04_plot_dreamlet_sn.R")

# slurmjobs::job_loop(
#     loops = list(model = names(dreamlet_models_sn)),
#     name = "04_plot_dreamlet_sn",
#     create_shell = TRUE,
#     create_script = FALSE
# )

