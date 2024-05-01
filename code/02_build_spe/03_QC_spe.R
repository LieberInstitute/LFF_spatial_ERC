
library("spatialLIBD")
library("tidyverse")
library("gridExtra")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "02_build_spe", "03_QC_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load spe
spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))

#### UMI ####
summary(spe$sum_umi[spe$inTissue])

## in tissue

## out of tissue

