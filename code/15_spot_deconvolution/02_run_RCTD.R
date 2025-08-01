## Louise Huuki-Myers, July 2025
## Run RCTD spot deconvolution

#### set up ####
library("SpatialExperiment")
library("SummarizedExperiment")
library("spacexr")
library("here")
library("sessioninfo")
library("getopt")

# Import command-line parameters
scec <- matrix(
    c("cell_type_col", "c", "1", "character", "cell type column",
      "sample", "s", "1", "character", "sample to subset"),
    ncol = 5, byrow = TRUE
)
opt <- getopt(scec)

#test
# opt$cell_type_col = "cell_type_broad"
# opt$sample = "Br5415"

print(opt)

# plot_dir <- here("plots", "15_spot_deconvolution", "01_runRCTD")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "15_spot_deconvolution", "02_run_RCTD", opt$cell_type_col)
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load RCTD data ####
rctd_data <- readRDS(here("processed-data", "15_spot_deconvolution", "01_prep_RCTD", sprintf("rctd_data_%s.rds", opt$cell_type_col)))

## subset to one sample
rctd_data$spatial_experiment <- rctd_data$spatial_experiment[,rctd_data$spatial_experiment$sample_id == opt$sample]
message(sprintf("Subset to sample_id = %s, ncol = %i", opt$sample, ncol(rctd_data$spatial_experiment)))

message("n cell types: ", length(rctd_data$cell_type_info$info[[2]]))

#### Run RTCD  ####
message(Sys.time(), " - Run RCTD")
## max_multi_types = median(spe$CNmask_dark_blue)
results_spe <- runRctd(rctd_data, rctd_mode = "multi", max_multi_types = 5, max_cores =8)

results_spe$sample_id

message(Sys.time(), " - Done.. saving results")
saveRDS(results_spe, file = here(data_dir, sprintf("RCTD_%s-%s.rds", opt$cell_type_col, opt$sample)))

# samples <- sort(unique(rctd_data$spatial_experiment$sample_id))

# slurmjobs::job_loop(loops = list(cell_type_col = c("cell_type_broad", "cell_type_anno"),
#                                  sample = samples),
#                     create_shell = TRUE,
#                     name = "02_run_RCTD",
#                     create_script = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
