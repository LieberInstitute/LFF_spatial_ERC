# Aug, 2024 - Louise Huuki-Myers
# Read in output from 01_get_droplet_scores, filter out empty droplets
# Filter for other quality metrics, build initial sce object

library("SingleCellExperiment")
library("DropletUtils")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")
library("EnsDb.Hsapiens.v86")
library("scran")
library("scuttle")
library("scater")
library("scDblFinder")
library("HDF5Array")

plot_dir <- here("plots", "04_snRNA-seq", "02_droplet_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "04_snRNA-seq", "02_droplet_QC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Sample names 
sample_list <- list.files(here("processed-data", "03_cellranger"))
length(sample_list)

#### get sample info ####
sample_info <- read_csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

#### get droplet scores ####
droplet_counts <- tibble(chromium_id = character(), n_pre_drop = numeric(), n_post_drop = numeric())
# droplet_counts <- data.frame(colnames = c("Sample", "pre_drop", "post_drop"))

message("Processing Droplet scores...")
sce_list <- map(sample_list, function(sample){
    fn10x <-  here("processed-data", "03_cellranger", sample, "outs", "raw_feature_bc_matrix") 
    fnDropScore <- here("processed-data", "04_snRNA-seq", "01_get_droplet_scores", paste0("droplet_scores_", sample, ".Rdata"))
    
    sce <- read10xCounts(fn10x, col.names=TRUE)
    load(fnDropScore)
    ncol_preDrop = ncol(sce)
    # 1148322
    
    sce <- sce[, which(e.out$FDR <= 0.001)]
    ncol_postDrop = ncol(sce)
    # 3622
    
    message(Sys.time(), " - ", sample, ": Pre-Drop = ", ncol_preDrop, " and Post-Drop = ", ncol_postDrop)
    ## save pre-post drops to table
    droplet_counts <<- droplet_counts |> add_row(chromium_id = sample, n_pre_drop = ncol_preDrop, n_post_drop = ncol_postDrop)
    
    return(sce)
})

## What were the knee cutoffs?
log_fn <- list.files(here("code", "04_snRNA-seq", "logs"), pattern = "01_get_droplet_scores_", full.names = TRUE)
logs <- map(log_fn, readLines)
names(logs) <- gsub("_[0-9]+.txt", "", gsub("01_get_droplet_scores_", "", basename(log_fn)))

knee_lower <- map_dbl(logs, ~ parse_number(.x[grepl("knee_lower =", .x)]))

# get manual knee
use_manual <- map_lgl(logs, ~ any(grepl("Use manual knee", .x)))
knee_manual <- as.integer(map2(logs, use_manual, ~ifelse(.y, parse_number(.x[grepl("Use manual knee", .x)]), NA)))

identical(names(knee_lower), droplet_counts$chromium_id)

droplet_counts$knee_lower <- knee_lower
droplet_counts$knee_manual <- knee_manual

summary(droplet_counts)
# chromium_id          n_pre_drop       n_post_drop     knee_lower     knee_manual 
# Length:31          Min.   : 581462   Min.   :2169   Min.   :206.0   Min.   :200  
# Class :character   1st Qu.:1048624   1st Qu.:4676   1st Qu.:223.5   1st Qu.:200  
# Mode  :character   Median :1262252   Median :5222   Median :273.0   Median :200  
#                    Mean   :1241313   Mean   :5185   Mean   :353.8   Mean   :200  
#                    3rd Qu.:1411650   3rd Qu.:5815   3rd Qu.:313.0   3rd Qu.:200  
#                    Max.   :1732806   Max.   :7847   Max.   :984.0   Max.   :200  
#                                                                     NA's   :25 

## total droplets sequenced
sum(droplet_counts$n_pre_drop)
# [1] 38480710

sum(droplet_counts$n_post_drop)
# [1] 160743

table(is.na(knee_manual))
# FALSE  TRUE 
# 6    25 

#### Compile sample quality info ####
message("Compile Sample Quality Info...")
## retrieve median reads per cell
cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger"))
cell_ranger_metrics <- map_dfr(cell_ranger_out_paths, 
                               ~read_csv(here("processed-data", "03_cellranger", .x, "outs", "metrics_summary.csv"),
                                         show_col_types = FALSE) |>
                                   mutate(chromium_id = .x, .before = 1))

summary(cell_ranger_metrics)
# chromium_id        Estimated Number of Cells Mean Reads per Cell Median Genes per Cell Number of Reads    
# Length:31          Min.   :1979              Min.   : 62782      Min.   :1416          Min.   :377788420  
# Class :character   1st Qu.:4756              1st Qu.: 84448      1st Qu.:2202          1st Qu.:465472465  
# Mode  :character   Median :5325              Median :101847      Median :2452          Median :525073680  
#                    Mean   :5268              Mean   :109874      Mean   :2515          Mean   :541737761  
#                    3rd Qu.:6069              3rd Qu.:120272      3rd Qu.:2847          3rd Qu.:578512992  
#                    Max.   :7830              Max.   :228840      Max.   :3613          Max.   :793042874  
# Valid Barcodes     Sequencing Saturation Q30 Bases in Barcode Q30 Bases in RNA Read Q30 Bases in UMI  
# Length:31          Length:31             Length:31            Length:31             Length:31         
# Class :character   Class :character      Class :character     Class :character      Class :character  
# Mode  :character   Mode  :character      Mode  :character     Mode  :character      Mode  :character  
# 
# 
# 
# Reads Mapped to Genome Reads Mapped Confidently to Genome Reads Mapped Confidently to Intergenic Regions
# Length:31              Length:31                          Length:31                                     
# Class :character       Class :character                   Class :character                              
# Mode  :character       Mode  :character                   Mode  :character                              
# 
# 
# 
# Reads Mapped Confidently to Intronic Regions Reads Mapped Confidently to Exonic Regions
# Length:31                                    Length:31                                 
# Class :character                             Class :character                          
# Mode  :character                             Mode  :character                          
# 
# 
# 
# Reads Mapped Confidently to Transcriptome Reads Mapped Antisense to Gene Fraction Reads in Cells Total Genes Detected
# Length:31                                 Length:31                      Length:31               Min.   :29063       
# Class :character                          Class :character               Class :character        1st Qu.:33028       
# Mode  :character                          Mode  :character               Mode  :character        Median :33730       
# Mean   :33428       
# 3rd Qu.:34162       
# Max.   :34831       
# Median UMI Counts per Cell
# Min.   : 2530             
# 1st Qu.: 4538             
# Median : 5369             
# Mean   : 5797             
# 3rd Qu.: 7141             
# Max.   :10275             

write_csv(cell_ranger_metrics, file = here(data_dir, "cell_ranger_metrics.csv"))

## create QC summary
sample_qc_summary <- sample_info |> 
    dplyr::select(BrNum, chromium_id, exp_round, seq_round, Ancestry, APOE, APOE_carrier) |>
    left_join(cell_ranger_metrics |> 
                  dplyr::select(1:6)) |> 
    left_join(droplet_counts) 

summary(sample_qc_summary)

#### Plot  drop values ####
message("Plot N nuclei after excluding empty droplets...")
# n nuclei barplot
droplet_barplot <- sample_qc_summary |>
    mutate(Sample = paste0(chromium_id, "-", BrNum)) |>
    ggplot(aes(x = reorder(Sample, n_post_drop), y = n_post_drop)) +
    geom_col() +
    geom_text(aes(label = n_post_drop), nudge_y = -290, color = "white") +
    theme_bw() +
    coord_flip() +
    labs(x = "Sample", y = "Non-empty Droplets")

ggsave(droplet_barplot, filename = here(plot_dir, "erc_n_droplets_post_drop_barplot.png"))

# boxplot by seq round
droplet_boxplot_seq <- sample_qc_summary |>
    ggplot(aes(x = as.factor(seq_round), y = n_post_drop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2) + 
    theme_bw() +
    labs(x = "Sequencing Round", y = "Non-empty Droplets")

ggsave(droplet_boxplot_seq, filename = here(plot_dir, "erc_n_droplets_post_drop_boxplot_seq.png"))

#boxplot by APOE
droplet_boxplot_apoe <- sample_qc_summary |>
    ggplot(aes(x = APOE, y = n_post_drop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2) + 
    theme_bw() +
    labs(y = "Non-empty Droplets")

ggsave(droplet_boxplot_apoe, filename = here(plot_dir, "erc_n_droplets_post_drop_boxplot_apoe.png"))

droplet_boxplot_apoe_carrier <- sample_qc_summary |>
    ggplot(aes(x = APOE_carrier, y = n_post_drop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2) + 
    theme_bw() +
    labs(y = "Non-empty Droplets")

ggsave(droplet_boxplot_apoe_carrier, filename = here(plot_dir, "erc_n_droplets_post_drop_boxplot_apoe_carrier.png"))

## compare droplet scores to cell ranger n cell predictions
post_drop_v_est_cells <- sample_qc_summary |>
    ggplot(aes(n_post_drop, `Estimated Number of Cells`)) +
    geom_point() +
    geom_text_repel(aes(label = chromium_id), size = 2) +
    geom_abline(color = "red") +
    labs(x = "n non-empty droplets (DropletUtils)", "n cells (CellRanger)")

ggsave(post_drop_v_est_cells, filename = here(plot_dir, "erc_post_drop_v_est_cells.png"))

## compare cells vs. mean reads
post_drop_v_mean_reads <- sample_qc_summary |>
    ggplot(aes(n_post_drop, `Mean Reads per Cell`)) +
    geom_point() +
    geom_text_repel(aes(label = chromium_id), size = 2) +
    labs(x = "n non-empty droplets (DropletUtils)") +
    geom_vline(xintercept = 9000, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 50000, color = "blue", linetype = "dashed")

ggsave(post_drop_v_mean_reads, filename = here(plot_dir, "erc_post_drop_v_mean_reads.png"))

#### Build SCE ####
message("Build Initial SCE object...")
table(map_int(sce_list, nrow))

## check row data is identical
all(map_lgl(sce_list, ~identical(rownames(sce_list[[1]]), rownames(.x))))

## rbind to one sce
sce <- do.call("cbind", sce_list)
dim(sce)
# [1]  38606 160743

#### Build colData ####
message("Build colData...")
colData(sce)

sce$path <- sce$Sample
sce$chromium_id <- gsub(".*?([0-9]+c_ERC_SVB)/outs.*?$", "\\1", sce$Sample)

pd <- as.data.frame(colData(sce)) |>
    left_join(sample_info) |>
    mutate(Barcode = paste0(Barcode, "_", BrNum),
           exp_round = factor(exp_round),
           seq_round = factor(seq_round),
           Ancestry = factor(Ancestry),
           APOE = factor(APOE),
           APOE_carrier = factor(APOE_carrier),
           Sex = factor(Sex)
           ) |>
    dplyr::select(Barcode, sample_id, chromium_id, BrNum, APOE:seq_round, path)

colData(sce) <- DataFrame(pd)

#### Compute QC metrics ####
message("\nCompute QC metrics...")

location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce)$ID,
                   column = "SEQNAME", keytype = "GENEID")
# Unable to map 4792 of 38606 requested IDs.
sce <- scuttle::addPerCellQC(
    sce,
    subsets = list(Mito = which(location=="MT"))
)

table(sce$sum < 200)

#### Doublet detection #### 
message("\nDoublet Detection...")
set.seed(821)
message(Sys.time(), " - ", "Run scDblFinder")
sce <- scDblFinder(sce, samples = sce$sample_id)
message(Sys.time(), " - Done")

table(sce$scDblFinder.class)
# singlet doublet 
# 149923   10820

100*sum(sce$scDblFinder.class == "doublet")/ncol(sce)
# 10X-like data tends to have roughly 1% per 1000 cells captured

#### Visualize doublet scores ####
dbl_violin <- plotColData(sce, x = "BrNum", y = "scDblFinder.score", colour_by = "scDblFinder.class") +
    facet_wrap(~ sce$seq_round, scales = "free_x", nrow = 1) +
    theme_bw()

ggsave(dbl_violin, filename = here(plot_dir, "erc_sn_doublet_scores_violin.png"), width = 21)

dbl_df <- colData(sce) |>
    as.data.frame() |>
    dplyr::select(sample_id, scDblFinder.score, scDblFinder.class)

doubletScore_summary <- dbl_df |>
    group_by(sample_id) |>
    dplyr::summarize(
        doubletScore_median = median(scDblFinder.score),
        n_doublets = sum(scDblFinder.class == "doublet")
    ) 

summary(doubletScore_summary)


#### find qc metric outliers w/ Scuttle::isOutlier ####
message("\nFind QC metric outliers...")

sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$sample_id)
sce$low_sum <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$sample_id)
sce$low_detected <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$sample_id)

## drop nuclei that are outliers in any of the three metrics
sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected
table(sce$discard_auto)
# FALSE   TRUE 
# 140119  20624

addmargins(table(sce$BrNum, sce$discard_auto))

## add by batch for ref
sce$high_mito_batch <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$seq_round)
sce$low_sum_batch <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$seq_round)
sce$low_detected_batch <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$seq_round)

#### QC plots ####
message("\nQC plots...")

qc_metrics <- c("sum", "detected", "subsets_Mito_percent")
names(qc_metrics) <- qc_metrics
qc_cutoff <- c("low_sum", "low_detected", "high_mito")

qc_violin_plots <- map2(qc_metrics, qc_cutoff,
                       ~plotColData(sce, x = "BrNum", y = .x, colour_by = .y) +
                           # scale_y_log10() +
                           # ggtitle(.x) +
                           facet_wrap(~ sce$seq_round, scales = "free_x", nrow = 1) +
                           theme_bw())

qc_violin_plots$sum <- qc_violin_plots$sum + scale_y_log10() 
qc_violin_plots$detected <- qc_violin_plots$detected + scale_y_log10()

walk2(qc_violin_plots, names(qc_violin_plots),
      ~ggsave(.x, filename = here(plot_dir, paste0("erc_sn_QC_outlier-",.y ,".png")),
              width = 21))

qc_detected_v_mito <- plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
) + theme_bw()

ggsave(qc_detected_v_mito, filename = here(plot_dir, "erc_sn_QC_detected_v_mito.png"))

# Detected features vs total count
qc_sum_v_detected <-plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
) + theme_bw()

ggsave(qc_sum_v_detected, filename = here(plot_dir, "erc_sn_QC_sum_v_detected.png"))

pdf(here(plot_dir, "erc_sn_QC_outliers.pdf"), width = 21)
print(qc_violin_plots)
print(qc_detected_v_mito)
print(qc_sum_v_detected)
dev.off()

#### Add qc counts to qc summary ####
message("\nAdd qc counts to qc summary...")

sample_qc_summary <- 
    sample_qc_summary |>
    left_join(doubletScore_summary |> rename(BrNum = sample_id)) |> ## add doublet data
    left_join(as.data.frame(colData(sce)) |> 
                  group_by(BrNum) |>
                  dplyr::summarise(n_low_sum = sum(low_sum),
                            n_low_detected = sum(low_detected), 
                            n_high_mito = sum(high_mito),
                            n_discard_auto = sum(discard_auto),
                            n_postQC = n()-n_discard_auto))

              

## summary of n post-QC
sum(sample_qc_summary$n_postQC) # [1] 140119
100*sum(sample_qc_summary$n_postQC)/ncol(sce) 
# [1] 87.16958
summary(sample_qc_summary$n_postQC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1811    4010    4589    4520    5100    6946 

sample_qc_summary |>
    group_by(APOE_carrier) |>
    summarise(total = sum(n_postQC),
              median = median(n_postQC))
# APOE_carrier total median
# <chr>        <int>  <dbl>
# 1 E2+          66605  4644.
# 2 E4+          73514  4250 

#### plot post QC values ####
message("\nPlot post QC values...")

n_nuc_qc_barplot <- sample_qc_summary |>
    mutate(BrNum = fct_reorder(BrNum, n_post_drop)) |>
    dplyr::select(BrNum, n_doublets, n_discard_auto, n_postQC) |>
    pivot_longer(!BrNum, values_to = "n_nuclei", names_to = "QC_class") |>
    group_by(BrNum) |>
    mutate(QC_class = factor(QC_class, levels = c("n_doublets", "n_discard_auto", "n_postQC"))) |>
    arrange(desc(QC_class)) |>
    mutate(text_y = cumsum(n_nuclei)) |>
    ggplot(aes(x = BrNum, y = n_nuclei, fill = QC_class)) +
    geom_col() +
    geom_text(aes(label = n_nuclei, y = text_y), nudge_y = -200, color = "white", size = 2) +
    theme_bw() +
    coord_flip() +
    labs(x = "Sample", y = "N Nuclei")

ggsave(n_nuc_qc_barplot, filename = here(plot_dir, "erc_n_nuc_post_QC_barplot.png"))


postQC_boxplot_seq <- sample_qc_summary |>
    ggplot(aes(x = seq_round, y = n_postQC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2) + 
    theme_bw() +
    labs(x = "Sequencing Round", y = "Nuclei Post-QC")

ggsave(postQC_boxplot_seq, filename = here(plot_dir, "erc_n_post_QC_boxplot_seq.png"))

postQC_boxplot_APOE_carrier <- sample_qc_summary |>
    ggplot(aes(x = APOE_carrier, y = n_postQC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2, color = "grey40") + 
    theme_bw() +
    labs(y = "Nuclei Post-QC")

ggsave(postQC_boxplot_APOE_carrier, filename = here(plot_dir, "erc_n_post_QC_boxplot_APOE_carrier.png"), height = 5, width = 5)

#### Drop Auto-drop nuclei ####
message("\nDrop Nuclei from sce...")
sce <- sce[,!sce$discard_auto]
dim(sce)
# [1]  38606 140119

table(sce$BrNum)

# check n per sample match
all(sample_qc_summary |> arrange(BrNum) |> pull(n_postQC) == table(sce$BrNum))

####  get logcounts ####
message(Sys.time(), " - logNormCounts")
sce <- logNormCounts(sce)

#### Save Data ####
message("\nSave Data...")
write_csv(sample_qc_summary, file = here("processed-data", "04_snRNA-seq", "02_droplet_QC", "erc_sn_sample_QC_summary.csv"))

## save HDF5
message(Sys.time(), " - Save hdf5")
saveHDF5SummarizedExperiment(sce,
                             dir = here("processed-data", "sce_objects", "hdf5_sce_postQC"), prefix = "", replace = TRUE,
                             chunkdim = NULL, level = NULL, as.sparse = TRUE,
                             verbose = TRUE)
message(Sys.time(), " - Done save hdf5")

# slurmjobs::job_single('02_droplet_QC', create_shell = TRUE, memory = '25G', command = "Rscript 02_droplet_QC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
