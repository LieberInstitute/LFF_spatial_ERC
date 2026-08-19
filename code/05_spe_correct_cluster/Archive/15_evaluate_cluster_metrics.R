## October 2024, Louise Huuki-Myers
## Evaluate fast h+ metric for spatial clusters

library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

data_dir <- here("processed-data", "05_spe_correct_cluster", "14_fasthplus")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load HDF5 SPE")
spe <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "spe_objects", "spe_ERC"))

#### fast H+ data ####

fasthplus <- read.csv(here("processed-data","05_spe_correct_cluster","14_fasthplus","fasthplus_results.csv")) |> 
    filter(!grepl("cluster",cluster.fasthplus)) |>
    separate(cluster.fasthplus, into = c(NA, "cluster", "fasthplus"), sep = " ") |>
    as_tibble() |>
    mutate(k = parse_number(cluster)) |>
    arrange(k)

fasthplus |> arrange(desc(fasthplus))

## get fasthplus from logs
fh_logs <- map(list.files(here("code", "05_spe_correct_cluster", "logs"), pattern = "14_fasthplus_BayesSpace_SVGm", full.names = TRUE),
               readLines)

fh_text <- map(fh_logs, ~.x[grep("Results: ", .x)])
parse_number("Results: BayesSpace_SVGm_k22H+ :0.228198132474115")




