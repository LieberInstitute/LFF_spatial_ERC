# Aug, 2024 - Louise Huuki-Myers
# Read in output from 01_get_droplet_scores, filter out empty droplets
# Filter for other quality metrics, build initial sce object

library("SingleCellExperiment")
library("DropletUtils")
library("tidyverse")
library("ggrepel")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "04_snRNA-seq", "02_droplet_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "04_snRNA-seq", "02_droplet_QC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Sample names 
sample_list <- list.files(here("processed-data", "03_cellranger"))
length(sample_list)

#### get sample info ####
sample_info <- read_csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

sample_info <- sample_info |> 
    mutate(sample_id = gsub("Hs_", "", Brain))

colnames(sample_info)

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
    
    message(Sys.time(), "-", sample, ": Pre-Drop = ", ncol_preDrop, " and Post-Drop = ", ncol_postDrop)
    ## save pre-post drops to table
    droplet_counts <<- droplet_counts |> add_row(Sample = sample, pre_drop = ncol_preDrop, post_drop = ncol_postDrop)
    
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

## retrieve median reads per cell
cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger"))
cell_ranger_metrics <- map_dfr(cell_ranger_out_paths, 
                               ~read_csv(here("processed-data", "03_cellranger", .x, "outs", "metrics_summary.csv")) |>
                                   mutate(chromium_id = .x, .before = 1))

cell_ranger_metrics <- cell_ranger_metrics |> left_join(droplet_counts |> rename(chromium_id = Sample))


summary(cell_ranger_metrics)

write_csv(cell_ranger_metrics, file = here(data_dir, "cell_ranger_metrics.csv"))

droplet_counts |>
    arrange(post_drop)

droplet_boxplot <- droplet_counts |>
    ggplot(aes(x = reorder(Sample, post_drop), y = post_drop)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    # geom_hline(yintercept = 5000, color = "red", linetype = "dashed") +
    labs(x = "Sample", y = "Non-empty Droplets")

ggsave(droplet_boxplot, filename = here(plot_dir, "erc_n_droplets_post_drop.png"))

post_drop_v_est_cells <- cell_ranger_metrics |>
    ggplot(aes(post_drop, `Estimated Number of Cells`)) +
    geom_point() +
    geom_text_repel(aes(label = chromium_id), size = 2) +
    geom_abline(color = "red") +
    labs(x = "n non-empty droplets (DropletUtils)", "n cells (CellRanger)")

ggsave(post_drop_v_est_cells, filename = here(plot_dir, "erc_post_drop_v_est_cells.png"))


post_drop_v_mean_reads <- cell_ranger_metrics |>
    ggplot(aes(post_drop, `Mean Reads per Cell`)) +
    geom_point() +
    geom_text_repel(aes(label = chromium_id), size = 2) +
    labs(x = "n non-empty droplets (DropletUtils)") +
    geom_vline(xintercept = 9000, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 50000, color = "blue", linetype = "dashed")

ggsave(post_drop_v_mean_reads, filename = here(plot_dir, "erc_post_drop_v_mean_reads.png"))

#### Build SCE ####
map_int(sce_list, nrow)

## check row data is identical
all(map_lgl(sce_list, ~identical(rownames(sce_list[[1]]), rownames(.x))))

## rbind to one sce
sce_erc_preQC <- do.call("cbind", sce_list)
dim(sce_erc_preQC)
# [1]  38606 154072

## Save preQC-sce
save(sce_erc_preQC, file = here("processed-data", "sce_objects", 
                               "sce_erc_preQC.Rdata"))
