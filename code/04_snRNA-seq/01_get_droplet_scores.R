
library("SingleCellExperiment")
library("DropletUtils")
# library("BiocParallel")
library("scuttle")
library("tidyverse")
library("here")
library("getopt")
library("sessioninfo")

## check dirs 
data_dir <- here("processed-data", "04_snRNA-seq", "01_get_droplet_scores")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "04_snRNA-seq", "01_get_droplet_scores")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## get sample i
# Import command-line parameters
spec <- matrix(
    c(
        c("sample_id", "manual_knee"),
        c("s", "k"),
        c("1","2"),
        c("character", "logical"),
        c("snRNA_seq sample", "use manual knee=200")
    ),
    ncol = 5
)
opt <- getopt(spec)

## Default knee param
if(is.null(opt$manual_knee)) opt$manual_knee <- FALSE

sample = opt$sample_id

message("for sample:", sample, "use knee=200: ", opt$manual_knee)

sample_path <- here("processed-data", "03_cellranger", sample, "outs", "raw_feature_bc_matrix") 
stopifnot(file.exists(sample_path))

message(Sys.time(), " - Reading data from ", sample_path)

#### Load & Subset raw data ####
sce <- read10xCounts(sample_path, col.names=TRUE) 

message("ncol:", ncol(sce))

#### Run barcodeRanks to find knee ####
message(Sys.time(), " - Running barcode ranks.")

bcRanks <- barcodeRanks(sce, fit.bounds = c(10, 1e3))

knee_lower <- metadata(bcRanks)$knee + 100
message(
  "Second knee point = ", metadata(bcRanks)$knee, "\n",
  "knee_lower = ", knee_lower
)

if(knee_lower < 200) warning("knee_lower < 200: may have nuclei w/ sums too low for doublet detection, consider manual_knee=TRUE")

#### Use Manual Knee if supplied ####
if(opt$manual_knee){
    message("Use manual knee = 200") # 200 based on minimal reads for dbl detection
    knee_lower <- 200
}

#### Run emptyDrops w/ knee + 100 ####
set.seed(100)
message(Sys.time(), " - Starting emptyDrops")
e.out <- DropletUtils::emptyDrops(
  sce,
  niters = 30000,
  lower = knee_lower
  # ,
  # BPPARAM = BiocParallel::MulticoreParam(4)
)
message(Sys.time(), " - Done: saving data")

save(e.out, file = here(data_dir, paste0("droplet_scores_", sample, ".Rdata")))

#### QC Plots ####
message("QC check")
FDR_cutoff <- 0.001
addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))

n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
message(n_cell_anno)

droplet_elbow_plot <- as.data.frame(bcRanks) %>%
  add_column(FDR = e.out$FDR) %>%
  ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
  annotate("text", x = 10,
           y = metadata(bcRanks)$knee, 
           label = paste("Second Knee:", metadata(bcRanks)$knee),
           vjust = -1, 
           color = "gray") +
  geom_hline(yintercept = metadata(bcRanks)$knee + 100, linetype = "dashed") +
  annotate("text", x = 100, 
           y = metadata(bcRanks)$knee + 100, 
           label = paste("Knee est 'lower':", metadata(bcRanks)$knee + 100), 
           vjust = -0.5) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Barcode Rank",
    y = "Total UMIs",
    title = paste("Sample", sample),
    subtitle = n_cell_anno,
    color = paste("FDR <", FDR_cutoff)
  ) +
    theme_bw() +
    theme(legend.position = "bottom",
          text = element_text(size = 15)) 

## Add Manual Knee if used
if(opt$manual_knee){
    droplet_elbow_plot <- droplet_elbow_plot + 
        geom_hline(yintercept = knee_lower, linetype = "dashed", color = "red") +
        annotate("text", x = 10, y = opt$knee, label = paste("Manual Knee:", knee_lower), vjust = -1, color = "red") 
}

ggsave(droplet_elbow_plot, filename = here(plot_dir, paste0("droplet_qc_", sample, ".png")))

## slurmjobs defined in 00_sn_data_check.R

# slurmjobs::job_loop(
#     loops = list(sample_id = c()),
#     name = "01_get_droplet_scores_manual_knee",
#     create_shell = TRUE
# )

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()
