

library(here)

metadata_sn <- read.csv(here("processed-data", "00_project_prep", "02_get_online_metadata","metadata_sn_plan.csv"))
metadata_sn

cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger"))
# [1] "10c_ERC_SVB" "11c_ERC_SVB" "12c_ERC_SVB" "13c_ERC_SVB" "14c_ERC_SVB" "1c_ERC_SVB"  "4c_ERC_SVB"  "5c_ERC_SVB" 
# [9] "6c_ERC_SVB"  "7c_ERC_SVB"  "8c_ERC_SVB"  "9c_ERC_SVB" 

