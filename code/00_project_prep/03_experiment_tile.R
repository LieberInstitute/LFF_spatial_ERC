
library("tidyverse")
library("googlesheets4")
library("sessioninfo")
library("here")

list.files(here("processed-data", "00_project_prep", "02_get_online_metadata"))

plot_dir <- here("plots", "00_project_prep", "03_experiment_tile")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

metadata_visium_list <- read_csv(file = here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_visium_list.csv"))
metadata_visium_plan <- read_csv(file = here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_visium_plan.csv")) 
metadata_sn_plan <- read_csv(file = here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_sn_plan.csv"))

metadata_sn_plan |> count(snRNA_complete)

metadata_plan <- metadata_visium_plan |>
    select(BrNum, APOE) |>
    

metadata_visium_plan |>
    ggplot(aes(x = BrNum, y = "Visium"), color = "black", fill = "red") +
    geom_tile() +
    facet_wrap(~APOE, nrow = 1, scales = "free_x")