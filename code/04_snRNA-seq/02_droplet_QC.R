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
library(scran)
library(scuttle)
library(scater)

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
droplet_counts <- tibble(chromium_id = character(), pre_drop = numeric(), post_drop = numeric())
# droplet_counts <- data.frame(colnames = c("Sample", "pre_drop", "post_drop"))

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
    droplet_counts <<- droplet_counts |> add_row(chromium_id = sample, pre_drop = ncol_preDrop, post_drop = ncol_postDrop)
    
    return(sce)
})


summary(droplet_counts)
# Sample             pre_drop         post_drop   
# Length:31          Min.   : 581462   Min.   :1794  
# Class :character   1st Qu.:1048624   1st Qu.:4420  
# Mode  :character   Median :1262252   Median :5072  
#                    Mean   :1241313   Mean   :4970  
#                    3rd Qu.:1411650   3rd Qu.:5720  
#                    Max.   :1732806   Max.   :7847 

#### Compile sample qualit info ####
## retrieve median reads per cell
cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger"))
cell_ranger_metrics <- map_dfr(cell_ranger_out_paths, 
                               ~read_csv(here("processed-data", "03_cellranger", .x, "outs", "metrics_summary.csv")) |>
                                   mutate(chromium_id = .x, .before = 1))

cell_ranger_metrics <- sample_info |> 
    select(BrNum, chromium_id, exp_round, seq_round, Ancestry, APOE) |>
    left_join(cell_ranger_metrics) |> 
    left_join(droplet_counts) 

summary(cell_ranger_metrics)

## maybe write later...?
write_csv(cell_ranger_metrics, file = here(data_dir, "cell_ranger_metrics.csv"))

#### Plot pot drop values ####
#barplot
droplet_barplot <- cell_ranger_metrics |>
    ggplot(aes(x = reorder(BrNum, post_drop), y = post_drop)) +
    geom_col() +
    geom_text(aes(label = post_drop), nudge_y = -275, color = "white") +
    theme_bw() +
    coord_flip() +
    labs(x = "Sample", y = "Non-empty Droplets")

ggsave(droplet_barplot, filename = here(plot_dir, "erc_n_droplets_post_drop_barplot.png"))

# boxplot by seq round
droplet_boxplot_seq <- cell_ranger_metrics |>
    ggplot(aes(x = as.factor(seq_round), y = post_drop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2) + 
    theme_bw() +
    labs(x = "Sequencing Round", y = "Non-empty Droplets")

ggsave(droplet_boxplot, filename = here(plot_dir, "erc_n_droplets_post_drop_boxplot_seq.png"))

#boxplot by APOE
droplet_boxplot_apoe <- cell_ranger_metrics |>
    ggplot(aes(x = APOE, y = post_drop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_text_repel(aes(label = BrNum), size = 2) + 
    theme_bw() +
    labs(y = "Non-empty Droplets")

ggsave(droplet_boxplot_apoe, filename = here(plot_dir, "erc_n_droplets_post_drop_boxplot_apoe.png"))

## compare droplet scores to cell ranger n cell predictions
post_drop_v_est_cells <- cell_ranger_metrics |>
    ggplot(aes(post_drop, `Estimated Number of Cells`)) +
    geom_point() +
    geom_text_repel(aes(label = chromium_id), size = 2) +
    geom_abline(color = "red") +
    labs(x = "n non-empty droplets (DropletUtils)", "n cells (CellRanger)")

ggsave(post_drop_v_est_cells, filename = here(plot_dir, "erc_post_drop_v_est_cells.png"))

## compare cells vs. mean reads
post_drop_v_mean_reads <- cell_ranger_metrics |>
    ggplot(aes(post_drop, `Mean Reads per Cell`)) +
    geom_point() +
    geom_text_repel(aes(label = chromium_id), size = 2) +
    labs(x = "n non-empty droplets (DropletUtils)") +
    geom_vline(xintercept = 9000, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 50000, color = "blue", linetype = "dashed")

ggsave(post_drop_v_mean_reads, filename = here(plot_dir, "erc_post_drop_v_mean_reads.png"))

#### Build SCE ####
table(map_int(sce_list, nrow))

## check row data is identical
all(map_lgl(sce_list, ~identical(rownames(sce_list[[1]]), rownames(.x))))

## rbind to one sce
sce <- do.call("cbind", sce_list)
dim(sce)
# [1]  38606 154072

#### Build colData ####
colData(sce)

sce$path <- sce$Sample
sce$chromiums_id <- gsub(".*?([0-9]+c_ERC_SVB)/outs.*?$", "\\1", sce$Sample)

pd <- as.data.frame(colData(sce)) |>
    left_join(sample_info) |>
    mutate(Barcode = paste0(Barcode, "_", BrNum),
           exp_round = factor(exp_round),
           seq_round = factor(seq_round),
           Ancestry = factor(Ancestry),
           APOE = factor(APOE),
           Sex = factor(Sex)
           ) |>
    select(Barcode, sample_id, chromium_id, BrNum, APOE:seq_round, path)

colData(sce) <- DataFrame(pd)


#### Compute QC metrics ####

location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce)$ID, 
                   column = "SEQNAME", keytype = "GENEID")
# Unable to map 4792 of 38606 requested IDs. 

sce <- scuttle::addPerCellQC(
    sce,
    # subsets = list(Mito = list(Mito = which(location=="MT")))
    subsets = list(Mito = which(location=="MT"))
    # ,
    # BPPARAM = BiocParallel::MulticoreParam(4)
)

## find outliers with Scran
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$sample_id)
sce$low_sum <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$sample_id)
sce$low_detected <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$sample_id)

## drop with any fails
sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected
table(sce$discard_auto)
# FALSE   TRUE 
# 125682  28390

addmargins(table(sce$BrNum, sce$discard_auto))

#### QC plots ####
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


walk2(qc_violin_plots, names(qc_violin_plots),
      ~ggsave(.x, filename = here(plot_dir, paste0("erc_sn_QC_outlier-",.y ,".png")),
              width = 21))

qc_detected_v_mito <- plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
)

ggsave(qc_detected_v_mito, filename = here(plot_dir, "erc_sn_QC_detected_v_mito.png"))

# Detected features vs total count
qc_sum_v_detected <-plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
)

ggsave(qc_sum_v_detected, filename = here(plot_dir, "erc_sn_QC_sum_v_detected.png"))

pdf(here(plot_dir, "erc_sn_QC_outliers.pdf"), width = 21)
print(qc_violin_plots)
print(qc_detected_v_mito)
print(qc_sum_v_detected)
dev.off()



## Save preQC-sce
save(sce, file = here("processed-data", "sce_objects", 
                               "sce.Rdata"))
