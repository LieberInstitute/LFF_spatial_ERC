## September 2024, Louise Huuki-Myers
## Check quality metrics and marker genes expression of SpatialDomains

library("spatialLIBD")
library("here")
library("sessioninfo")

#### define dirs ####
data_dir <- here("processed-data", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_spe_correct_cluster", "11_SpD_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


#### specify and load data ####

data


