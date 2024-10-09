## October 2024, Louise Huuki-Myers
## Load sample-wise nnSVG data, average the ranks of the SVGs

library("SpatialExperiment")
library("HDF5Array")
library(tidyverse)
library("here")
library("sessioninfo")

#### load spe data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "13_gather_nnSVG")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


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

nnSVG_avg <- nnSVG_data |>
    group_by(gene_id, gene_name) |>
    summarize(nnsvg_avg_rank = mean(rank),
              n = n()) |>
    arrange(nnsvg_avg_rank)

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

nnSVG_avg |>
    ggplot(aes(x = n, y = nnsvg_avg_rank)) +
    geom_point()

# A tibble: 5,879 × 4
# Groups:   gene_id [5,879]

## genes we're rareley evalauted in all samples
nnSVG_avg |> ungroup() |> count(n == 31)
# `n == 31`     n
# <lgl>     <int>
# 1 FALSE      5614
# 2 TRUE        265

# slurmjobs::job_single('13_gather_nnSVG', create_shell = TRUE, memory = '10G', command = "Rscript 13_gather_nnSVG.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()s

