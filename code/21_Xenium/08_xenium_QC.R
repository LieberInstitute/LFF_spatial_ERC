## Louise Huuki-Myers, April 2026
## Build SummarizedExperiment Objetct for Xenium data
## Adapted from https://github.com/LieberInstitute/xenium_NAC/blob/54a6bd78ef1127d0e41f037bd6a080579dad9020/code/03_QC/01_qc.R

#### Set up ####

library("here")
library("SpatialExperiment")
library("qs2")
library("scater")
library("tidyverse")

data_dir <- here("processed-data", "21_Xenium", "08_xenium_QC")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "21_Xenium", "08_xenium_QC")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
message(Sys.time(), "- Loda data")
spe <- qs_read(here("processed-data", "21_Xenium", "07_xenium_build_spe","spe_xenium.qs2"))

spe
# dim: 541 473611

#### Remove empty cells ####
message("Remove empty cells:", sum(colSums(counts(spe)) == 0)) #659
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

spe$sum_gex_low <- isOutlier(spe$sum_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=3, type="lower")
sum(spe$sum_gex_low) #1647

spe$detected_gex_low  <- isOutlier(spe$detected_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=4, type="lower") 
sum(spe$detected_gex_low) # 2204

# 5) Cell area outliers (after removing GEX-empties)
spe$cell_area_outlier <- isOutlier(spe$cell_area, log=FALSE, batch=as.factor(spe$Sample), nmads=3, type="both")  
sum(spe$cell_area_outlier) #402

# global outliers flag
spe$global_outliers <- as.logical(spe$exclude_any_neg) |
    as.logical(spe$detected_gex_low) | 
    as.logical(spe$sum_gex_low) |
    as.logical(spe$cell_area_outlier)

table(spe$global_outliers)

spe$APOE_ancestry <- paste0(gsub("\\+", "_", spe$APOE_carrier), spe$Ancestry)

QC_summary <- as.data.frame(colData(spe)) |>
    group_by(Run, chip, BrNum, APOE_carrier, Ancestry, APOE_ancestry) |>
    summarize(
        nCells = n(),
        CellsRetained = sum(!global_outliers),
        Outliers = sum(global_outliers),
        PercentRemoved  = Outliers/(CellsRetained + Outliers)*100,
        exclude_any_neg = sum(exclude_any_neg),
        detected_gex_4MAD  = sum(detected_gex_4MAD),
        sum_gex_4MAD = sum(sum_gex_4MAD),
        cell_area_4MAD  = sum(cell_area_4MAD)
    ) 

#### QC plots ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

pd <- as.data.frame(colData(spe))

## n cells bar plot 
n_cell_barplot <- QC_summary |>
    ggplot(aes(x = BrNum, y = nCells, fill = APOE_ancestry)) +
    geom_col() +
    facet_wrap(~Run+chip, scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = APOE_carrier_colors)

ggsave(n_cell_barplot, filename = here(plot_dir, "xenium_n_cell_barplot.png"), height = 5)

n_cell_boxplot <- QC_summary |>
    ggplot(aes(x = chip, y = nCells)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = APOE_ancestry), width = 0.2) +
    facet_wrap(~Run, scales = "free_x", nrow = 1) +
    theme_bw()  +
    scale_color_manual(values = APOE_carrier_colors)

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


