## October 2024, Louise Huuki-Myers
## Load sample-wise nnSVG data, average the ranks of the SVGs

library("SpatialExperiment")
library("HDF5Array")
library("tidyverse")
library("here")
library("sessioninfo")

#### load spe data ####
message(Sys.time(), " - Load HDF5 SPE")
# spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "13_gather_nnSVG")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "13_gather_nnSVG")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load nnSVG outputs ####

nnSVG_files <- list.files(here("processed-data", "05_spe_correct_cluster", "12_nnSVG"), full.names = TRUE)
names(nnSVG_files) <- gsub("nnSVG_(Br[0-9]+).csv", "\\1", basename(nnSVG_files))

nnSVG_data <- do.call("rbind", map2(nnSVG_files, names(nnSVG_files), ~read.csv(.x, row.names = 1) |> mutate(sample_id = .y)))

head(nnSVG_data)

nnSVG_data |> count(sample_id) |> arrange(n)
#    sample_id    n
# 1     Br5276  304
# 2     Br5529  848
# 3     Br6263 1159
# 4     Br6321 1555
# 5     Br6161 1602

dlpfc_layers_top100 <-  read_csv(here("processed-data", "00_project_prep", "06_marker_genes", "dlpfc_layers_top100.csv"))
length(unique(dlpfc_layers_top100$ensembl)) # [1] 1129

## other gene sets
load(here("processed-data", "02_build_spe","01_preprocess_spe", "top_hvgs.Rdata"), verbose = TRUE)
load(here("processed-data", "05_spe_correct_cluster", "03_GLM_Harmony", "top_hdgs.rdata"), verbose = TRUE)

# hdg_vs_hvg_overlap <- intersect(top.hvgs$p1, top.hdgs$`2k`)
# length(hdg_vs_hvg_overlap) #[1] 405
# 
# SVG_vs_hvg_overlap <- intersect(top.hvgs$p1, top.svg)
# SVG_vs_hdg_overlap <- intersect(top.hdgs$`1k`, top.svg)
# 
# length(SVG_vs_hvg_overlap) #[1] 176
# length(SVG_vs_hdg_overlap) #[1] 464

nnSVG_avg <- nnSVG_data |>
    group_by(gene_id, gene_name) |>
    summarize(nnsvg_avg_rank = mean(rank),
              n = n()) |>
    arrange(nnsvg_avg_rank) |>
    ungroup() |>
    mutate(top.svg = (n > 5 & nnsvg_avg_rank < 1500),
           DLPFC_layer_marker = gene_id %in% dlpfc_layers_top100$ensembl,
           top.hvg = gene_id %in% top.hvgs$p1,
           top.hdg = gene_id %in% top.hdgs$`1k`) 

nnSVG_avg |> count(top.svg)

nnSVG_avg |> count(top.svg, DLPFC_layer_marker)

nnSVG_avg |> count(top.svg, DLPFC_layer_marker)
nnSVG_avg |> count(top.svg, top.hvg, top.hdg, DLPFC_layer_marker)

write_csv(nnSVG_avg, file = here(data_dir, "nnSVG_avg.csv"))

# gene_id         gene_name nnsvg_avg_rank     n
# <chr>           <chr>              <dbl> <int>
# 1 ENSG00000197971 MBP                 2       31
# 2 ENSG00000131095 GFAP                9.29    31
# 3 ENSG00000106809 OGN                10        1
# 4 ENSG00000123560 PLP1               10.5     31
# 5 ENSG00000118271 TTR                12.7     18
# 6 ENSG00000132639 SNAP25             19.6     31
# 7 ENSG00000154146 NRGN               26.5     31
# 8 ENSG00000142173 COL6A2             29        1
# 9 ENSG00000165507 DEPP1              29.5      2
# 10 ENSG00000171617 ENC1               29.6     31

nnSVG_n_vs_avg_rank <- nnSVG_avg |>
    ggplot(aes(x = n, y = nnsvg_avg_rank, color = top.svg)) +
    geom_point()

ggsave(nnSVG_n_vs_avg_rank, filename = here(plot_dir, "nnSVG_n_vs_avg_rank.png"))

nnSVG_n_vs_avg_rank_DLPFC <- nnSVG_avg |>
    ggplot(aes(x = n, y = nnsvg_avg_rank, color = DLPFC_layer_marker)) +
    geom_point()

ggsave(nnSVG_n_vs_avg_rank_DLPFC, filename = here(plot_dir, "nnSVG_n_vs_avg_rank-DLPFC_markers.png"))

nnSVG_n_histogram <- nnSVG_avg |>
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 1)

ggsave(nnSVG_n_histogram, filename = here(plot_dir, "nnSVG_n_histogram.png"))

## genes were rarely evaluated in all samples
nnSVG_avg |> ungroup() |> count(n == 31)
# `n == 31`     n
# <lgl>     <int>
# 1 FALSE      5614
# 2 TRUE        265

nnSVG_avg |> ungroup() |> count(n >= 20)


top.svg <- nnSVG_avg |> 
    ungroup() |>
    filter(n > 5,
           nnsvg_avg_rank < 1500) |>
    pull(gene_id)

length(top.svg)
# [1] 1454

save(top.svg, file = here(data_dir, "top_svg.Rdata"))


#### Load SPE data ####
message(Sys.time(), " - Loading SPE data")
spe <- loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_postQC"))

 
# slurmjobs::job_single('13_gather_nnSVG', create_shell = TRUE, memory = '10G', command = "Rscript 13_gather_nnSVG.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()s

