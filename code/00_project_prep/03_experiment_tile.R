
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
    mutate(assay = "Visium", complete = TRUE) |>
    bind_rows(metadata_sn_plan |>
                  select(BrNum, APOE, complete = snRNA_complete) |>
                  mutate(assay = "snRNA-seq"))

metadata_plan |> count(assay, complete)
# assay     complete     n
# <chr>     <lgl>    <int>
# 1 Visium    TRUE        32
# 2 snRNA-seq FALSE       17
# 3 snRNA-seq TRUE        14

experiment_tile <- metadata_plan |>
    ggplot(aes(x = BrNum, y = assay, fill = complete)) +
    geom_tile(color = "black") +
    facet_grid(.~APOE, scales = "free_x", space='free')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
ggsave(experiment_tile, filename = here(plot_dir, "experiment_tile.png"), height = 4, width = 9)

# slurmjobs::job_single('03_experiment_tile', create_shell = TRUE, memory = '5G', command = "Rscript 03_experiment_tile.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

