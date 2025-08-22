## Louise Huuki-Myers, Aug 2025
## Find overlap in DEG and LR data

#### Set up ####
library("tidyverse")
library("here")
library("sessioninfo")
library("getopt")
library("data.table")

# Import command-line parameters
scec <- matrix(
    c("datatype", "d", "1", "character", "Data type"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

## test
# opt$datatype = "sn_broad"
# opt$datatype = "sn_fine"
# opt$datatype = "Visium"

data_dir <- here("processed-data", "17_cell_cell_comm", "10_LR_DEG_check")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "17_cell_cell_comm", "10_LR_DEG_check")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# load(here("processed-data", "project_colors.Rdata"))

#### load DE data ####
DE_data_fn <- here("processed-data", "13_compile_DGE", "01_compile_DGE", opt$datatype, sprintf("DGE_results_carrier_%s.Rds", opt$datatype))
file.exists(DE_data_fn)


#### load LR data ####
liana_fn <- here("processed-data", "17_cell_cell_comm", "liana", "ranked_results", sprintf("%s.csv.gz", tolower(gsub("sn_", "", opt$datatype))))
file.exists(liana_fn)

liana_data <- fread(liana_fn)


