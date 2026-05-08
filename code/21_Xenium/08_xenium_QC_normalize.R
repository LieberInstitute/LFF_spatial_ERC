## Louise Huuki-Myers, April 2026
## Build SummarizedExperiment Objetct for Xenium data
## Adapted from https://github.com/LieberInstitute/xenium_NAC/blob/54a6bd78ef1127d0e41f037bd6a080579dad9020/code/03_QC/01_qc.R

#### Set up ####

library("here")
library("SpatialExperiment")
library("qs2")
library("scater")
library("tidyverse")

data_dir <- here("processed-data", "21_Xenium", "08_xenium_QC_normalize")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "08_xenium_QC_normalize")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## SET UP MULTICORE PARAM
ncores <- 4  # match Sys.getenv('SLURM_CPUS_ON_NODE')
bp <- MulticoreParam(workers = ncores) #or bp <- MulticoreParam(4)

#### load data ####
message(Sys.time(), "- Load xenium data")
spe <- qs_read(here("processed-data", "21_Xenium", "07_xenium_build_spe","spe_xenium.qs2"))

spe
# dim: 541 473611

colnames(spe) <- paste0(spe$BrNum, "_", spe$Barcode)
any(duplicated(colnames(spe)))

spe$APOE_ancestry <- paste0(gsub("\\+", "_", spe$APOE_carrier), spe$Ancestry)

sample_nCells_preQC <- as.data.frame(colData(spe)) |>
    group_by(BrNum, APOE_carrier, Ancestry, APOE_ancestry, Run, chip) |>
    summarize(nCells_preQC = n()) 

#### Remove empty cells ####
message("Remove empty cells:", sum(colSums(counts(spe)) == 0)) #659
table(spe$BrNum, colSums(counts(spe)) == 0)

spe <- spe[, colSums(counts(spe)) > 0]           

#### Find Negative control genes for scuttle subsets ####

#Negative control probes are sequences that should not bind anything and therefore measure false positive rates and specificity
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
sum(is_neg) #20
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
sum(is_neg2) #41
#Unassigned probes help measure noise and off-target activity because they do not match any gene codewords
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
sum(is_unassigned) #61

is_anyneg <- is_neg | is_neg2 | is_unassigned
table(is_anyneg)
# FALSE  TRUE 
# 419   122

is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel
table(is_GEX)
# FALSE  TRUE 
# 175   366 

table(is_GEX, is_anyneg)
#          is_anyneg
# is_GEX  FALSE TRUE
# FALSE    53  122
# TRUE    366    0

#### Evaluate QC metrics scuttle::isOuliter() ####
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        any_neg = is_anyneg,
                                                        GEX = is_GEX))
## check if any cells have all negative control counts
any(spe$subsets_any_neg_percent == 100)

# 1) Subset to GEX rows (panel genes)
spe <- spe[is_GEX, ]
message("subset to GEX rows: ", nrow(spe))

# 2) Remove empty cells (based on GEX-only counts)
spe$sum_gex <- colSums(counts(spe))
message("sum gex")
summary(spe$sum_gex)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1     267     468     613     779    5570 

spe$detected_gex <- colSums(counts(spe) > 0)
message("detected gex")
summary(spe$detected_gex)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   59.00   86.00   90.08  117.00  244.00

table(spe$detected_gex == 1)
# FALSE   TRUE 
# 472692    260 

# 3) Noise fraction from FULL features (already computed pre-subset)
spe$exclude_any_neg <- spe$subsets_any_neg_percent >= 25              # 324 cells


# 4) Adaptive total counts and detected genes outliers (GEX-only metrics), per sample

spe$sum_gex_low <- isOutlier(spe$sum_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=4, type="lower")
message("sum_gex_low: ", sum(spe$sum_gex_low)) #4161

spe$detected_gex_low  <- isOutlier(spe$detected_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=4, type="lower") 
message("detected_gex_low: ", sum(spe$detected_gex_low)) # 2204

# 5) Cell area outliers (after removing GEX-empties)
spe$cell_area_outlier <- isOutlier(spe$cell_area, log=FALSE, batch=as.factor(spe$Sample), nmads=3, type="both")  
message("cell_area_outlier: ", sum(spe$cell_area_outlier)) #402

# global outliers flag
spe$global_outliers <- as.logical(spe$exclude_any_neg) |
    as.logical(spe$detected_gex_low) | 
    as.logical(spe$sum_gex_low) |
    as.logical(spe$cell_area_outlier)

message("All global outliers")
table(spe$global_outliers)

## summarize cutoffs
QC_cutoff_summary <- map2(
    c("sum_gex_low", "detected_gex_low", "cell_area_outlier"),
    c("sum_gex", "detected_gex", "cell_area"),
    ~{
        metric <- .y       
        outlier_col <- .x  
        
        as.data.frame(colData(spe)) |>
            group_by(BrNum) |>
            filter(!(!!sym(outlier_col))) |>
            summarise(
                min = min(!!sym(metric)),
                max = max(!!sym(metric))
            ) |>
            rename_with(~paste0(metric, "_", .x), c(min, max))
    }
) |>
    reduce(left_join, by = "BrNum")


QC_cutoff_summary |> arrange(cell_area_max)


QC_cutoff_summary_global <- map2_dfr(
    c("sum_gex_low", "detected_gex_low", "cell_area_outlier"),
    c("sum_gex", "detected_gex", "cell_area"),
    ~{
        metric <- .y       
        outlier_col <- .x  
        
        as.data.frame(colData(spe)) |>
            # group_by(BrNum) |>
            filter(!(!!sym(outlier_col)), !global_outliers) |>
            summarise(
                metric = metric,
                min = min(!!sym(metric)),
                max = max(!!sym(metric))
            ) 
    }
)

#        metric       min       max
# 1      sum_gex 13.000000 4835.0000
# 2 detected_gex  7.000000  237.0000
# 3    cell_area  2.935156  401.7552

## add additonal stats to QC summary
QC_summary <- sample_nCells_preQC |> 
    left_join(as.data.frame(colData(spe)) |>
                  group_by(BrNum) |>
                  summarize(
                      nCells = n(),
                      CellsRetained = sum(!global_outliers),
                      Outliers = sum(global_outliers),
                      PercentRemoved  = Outliers/(CellsRetained + Outliers)*100,
                      exclude_any_neg = sum(exclude_any_neg),
                      sum_gex_low  = sum(sum_gex_low),
                      detected_gex_low = sum(detected_gex_low),
                      cell_area_outlier  = sum(cell_area_outlier)
                  ) ) |>
    left_join(QC_cutoff_summary |> select(BrNum, sum_gex_min, detected_gex_min, cell_area_min, cell_area_max))


summary(QC_summary$sum_gex_min)
summary(QC_summary$detected_gex_min)

write.csv(QC_summary, file = here(data_dir, "ERC_Xenium_QC_summary.csv"))



#### QC plots ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

pd <- as.data.frame(colData(spe))

## n cells bar plot 
n_cell_barplot <- QC_summary |>
    ungroup() |>
    select(BrNum, APOE_ancestry, Run, chip, CellsRetained, Outliers) |>
    pivot_longer(!c(BrNum, APOE_ancestry, Run, chip), names_to = "QC", values_to = "nCells") |>
    mutate(QC = fct_rev(QC)) |>
    ggplot(aes(x = BrNum, y = nCells, fill = APOE_ancestry, color = QC)) +
    geom_col() +
    facet_wrap(~Run+chip, scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = APOE_carrier_anc_colors) +
    scale_color_manual(values = c(Outliers = "red", CellsRetained = "black"))

ggsave(n_cell_barplot, filename = here(plot_dir, "xenium_n_cell_barplot.png"), height = 5)

n_cell_boxplot <- QC_summary |>
    ggplot(aes(x = chip, y = CellsRetained)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = APOE_ancestry), width = 0.2) +
    facet_wrap(~Run, scales = "free_x", nrow = 1) +
    theme_bw()  +
    scale_color_manual(values = APOE_carrier_anc_colors)

ggsave(n_cell_boxplot, filename = here(plot_dir, "xenium_n_cell_boxplot.png"), height = 5)


## QC metrics by sample 
qc_metrics <- c("sum_gex", "detected_gex", "cell_area")
names(qc_metrics) <- qc_metrics
qc_cutoff <- c("sum_gex_low", "detected_gex_low", "cell_area_outlier")

qc_violin_plots <- map2(qc_metrics, qc_cutoff,
                        ~plotColData(spe, x = "BrNum", y = .x, colour_by = .y) +
                            # scale_y_log10() +
                            # ggtitle(.x) +
                            facet_wrap(~ spe$chip, scales = "free_x", nrow = 1) +
                            theme_bw())

walk2(qc_violin_plots, names(qc_violin_plots),
      ~ggsave(.x, filename = here(plot_dir, paste0("erc_xenium_QC_outlier-",.y ,".png")),
              width = 21))


#### GGpairs plots ####
library(GGally)

gg_QC_plot <- ggpairs(pd, columns = c("sum_gex", "detected_gex", "cell_area"), aes(colour = chip, size = 1, alpha = 0.5)) + theme_bw()
ggsave(gg_QC_plot, filename = here(plot_dir, "xenium_QC_metrics_ggpairs_plot.png"), height = 12, width = 12)

gg_QC_plot_out <- ggpairs(pd, columns = c("sum_gex", "detected_gex", "cell_area"), aes(colour = global_outliers)) + theme_bw()
ggsave(gg_QC_plot_out, filename = here(plot_dir, "xenium_QC_metrics_ggpairs_plot_outliers.png"), height = 12, width = 12)

#### Spot Sweeper ####

walk2(c("sum_gex", "detected_gex", "cell_area"),
       c("sum_gex_low", "detected_gex_low", "cell_area_outlier"),
       ~SpotSweeper::plotQCpdf(
           spe,
           sample_id  = "BrNum",
           metric     = .x, 
           outliers   = .y,
           point_size = 0.5,
           stroke     = 0.5,
           # colors     = viridis(256), 
           fname      = here(plot_dir,sprintf("spatial_QC_%s.pdf", .x))
       )
)

#### echer ####
library(escheR)

escher_QC <- map(sort(unique(spe$BrNum)), 
                   ~make_escheR(spe[,spe$BrNum==.x]) |>
                       add_fill(var = "sum_gex",point_size=0.6) |>
                       add_ground("global_outliers",point_size=0.6) +
                       scale_color_manual(values = c(`TRUE` = "red", `FALSE` = NA))
)

ggsave(patchwork::wrap_plots(escher_QC, ncol = 4), filename = here(plot_dir, "escher_QC.png"), height = 16, width = 16)



#### Vis_clus ####
# library(spatialLIBD)

head(spatialCoords(spe))
#      x_centroid y_centroid
# [1,]   642.6926   1825.158
# [2,]   584.5524   1850.862
# [3,]   661.8859   1807.252

# Error: Abnormal spatial coordinates: should have 'pxl_row_in_fullres' and 'pxl_col_in_fullres' columns.

# imgData(spe)
# 
# test <- vis_clus(spe,
#                  clustervar = "global_outliers",
#                  is_stitched = TRUE,
#                  spatial = FALSE)
# 
# test <- vis_image(spe)
# 
# 
# SpatialExperiment::scaleFactors(spe, sample = "sample01.1")

#### Drop Global Outliers ####

#Remove all of the low quality cells
table(spe$global_outliers)
spe <- spe[,!spe$global_outliers]
dim(spe)

#### Normalize ####
## adapted from https://github.com/LieberInstitute/xenium_NAC/blob/54a6bd78ef1127d0e41f037bd6a080579dad9020/code/04_normalization/01_area_based_norm.R
## Goal: Calculate size factors based on nucleus and cell area

pd <- as.data.frame(colData(spe))

#Based on Atta et al. 2024, generate non-count based scaling factors
#Calculate both cell area and nucleus area values
message("cell_area.sf")
spe$cell_area.sf <- spe$cell_area / median(spe$cell_area)
summary(spe$cell_area.sf)
message("nucleus_area.sf")
spe$nucleus_area.sf <- spe$nucleus_area / median(spe$nucleus_area)
summary(spe$nucleus_area.sf)


#Generate histograms 
cell_area_hist <- ggplot(colData(spe),aes(x = cell_area.sf)) +
    geom_histogram(bins = 50,fill = "#70bbfe",color = "#FFF") +
    labs(x = "Scaling Factor",
         y = "Frequency",
         title = "Cell area") +
    theme_bw() +
    theme(plot.title = element_text(hjust= 0.5))

ggsave(plot = cell_area_hist, filename = here(plot_dir,"cell_area_scaling_hist.pdf"))

nuc_area_hist <- ggplot(colData(spe),aes(x = nucleus_area.sf)) +
    geom_histogram(bins = 50,fill = "#70bbfe",color = "#FFF") +
    labs(x = "Scaling Factor",
         y = "Frequency",
         title = "Nucleus area") +
    theme_bw() +
    theme(plot.title = element_text(hjust= 0.5))

ggsave(plot = nuc_area_hist, filename = here(plot_dir,"nucleus_area_scaling_hist.pdf"))

# From https://github.com/LieberInstitute/spatialAmygdala/blob/77670e73360a112945a4b8257557194f9212bdb6/code/Xenium/03_quality_control/perCellQC/01_perCellQC.R#L102-L103
# normalize the counts by the nucleus and cell area scaling factors
message(Sys.time(), " -  Normalize by size factors")
assay(spe, "logcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$nucleus_area.sf, transform="log", assay.type="counts") ## segmentation of nuclei is more relable
assay(spe, "cell_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$cell_area.sf, transform="log", assay.type="counts")

#### Reduced dims ####

# use logcounts of your preferred normalization
# spe <- logNormCounts(spe, assay.type = "counts")

assay(spe, "logcounts") <- log2(assay(spe, "nucleus_normcounts") + 1)

message(Sys.time(), " - running PCA")
spe <- runPCA(spe,
              ncomponents = 30,
              BSPARAM = BiocSingular::IrlbaParam(),
              BPPARAM = bp
)

message(Sys.time(), " - running TSNE")
spe <- runTSNE(spe, dimred = "PCA", BPPARAM = bp)

message(Sys.time(), " - running UMAP")
spe <- runUMAP(spe, dimred = "PCA", BPPARAM = bp)

# Print reduced dimension names
message("\n\n", "Reduced Dim Names:\n")
reducedDimNames(spe)


#Save cleaned SPE
message(Sys.time(), " - Saving cleaned SPE object")
qs2::qs_save(spe, here(data_dir, "spe_xenium_QC.qs2"))

# elbow plot
var_explained <- data.frame(PC = seq_along(pca_var_pct),
                            var_explained = pca_var_pct) 

var_explained_elbow <- var_explained |>
    ggplot(aes(x = PC, y = var_explained)) +
    geom_point() +
    geom_line() +
    labs(x = "PC", y = "Variance Explained (%)") +
    theme_bw()

ggsave(var_explained_elbow, filename = here(plot_dir, "var_explained_elbow.png"))

source(here("code", "utils", "my_plot_reduced_dim.R"))

#### plot ####
## categorical
walk(c("sample_id", "seq_round", "chip","APOE"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "cat"))
walk(c("sample_id", "seq_round", "chip","APOE"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "TSNE", my_var = .x, var_type = "cat"))

## continuous
walk(c("sum_gex", "detected_gex", "cell_area"), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", my_var = .x, var_type = "con"))

my_genes <- c(Astro = "AQP4", 
              # Micro = "CD86",
              Micro = "CTSH",
              Oligo = "MBP",
              Oligo.M = "OPALIN",
              Oligo.3 = "LINGO2",
              OPC = "PDGFRA",
              Vasc = "PECAM1",
              Excit = "SLC17A7",
              Inhib = "GAD1")

walk2(my_genes, names(my_genes), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", var_type = "express", my_var = .x, suffix = .y, NA_gray = TRUE))
walk2(my_genes, names(my_genes), ~my_plot_reduced_dim(spe, prefix = "ERC_xenium", dimred = "UMAP", var_type = "express", my_var = .x, suffix = .y, NA_gray = TRUE))


# spe <- qs_read(here("processed-data", "21_Xenium", "08_xenium_QC_normalize","spe_xenium_QC.qs2"))


# slurmjobs::job_single('08_xenium_QC_normalize', create_shell = TRUE, memory = '100G', command = "Rscript 08_xenium_QC_normalize.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


