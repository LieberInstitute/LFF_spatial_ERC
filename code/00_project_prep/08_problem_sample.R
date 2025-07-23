## Louise Huuki-Myers, July 2025
## organize data and files relevant to problem donor Br1289

library(here)
library(tidyverse)

## read in sample info
sn_sample_info <- read.csv(here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

visium_sample_info <- read.csv(here("processed-data", "02_build_spe", "sample_info.csv"))


sn_id <- sn_sample_info |> filter(BrNum == "Br1289") |> pull(chromium_id)
# "19c_ERC_SVB"

visium_id <- visium_sample_info |> filter(BrNum == "Br1289") |> pull(VNum)
# "V13B23-364_D1"

sn_bamfile <- here("processed-data", "03_cellranger", sn_id, "outs", "possorted_genome_bam.bam")
# "/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/03_cellranger/19c_ERC_SVB/outs/possorted_genome_bam.bam"

visium_bamfile <- here("processed-data", "01_spaceranger", visium_id, "outs", "possorted_genome_bam.bam")
# "/dcs05/lieber/marmaypag/LFF_spatialERC_LIBD4140/LFF_spatial_ERC/processed-data/01_spaceranger/V13B23-364_D1/outs/possorted_genome_bam.bam"