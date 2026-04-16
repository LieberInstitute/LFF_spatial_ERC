## Louise Huuki-Myers, April 2026
## Check polygenic risk score for dononrs

library("tidyverse")
library("here")

data_dir <- here("processed-data", "00_project_prep", "12_polygenic_risk_score")
if(!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("processed-data", "00_project_prep", "12_polygenic_risk_score")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

