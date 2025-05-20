## Louise Huuki-Myers, May 2025
## Read and check ancestry data for donors

library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "00_project_prep", "04_ancestry_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## Feb 2025 update
## Shizhong FLARE estimates
ancestry_tab <- read.csv("/dcs04/lieber/shared/statsgen/local_ancestry/flare_array_data/global_ancestry.csv", row.names = 1) |>
    select(BrNum = ID, Afr = YRI, Eur = CEU)

head(ancestry_tab)

## Old data 
# ancestry_tab <- read.table(here("raw-data", "ancestry", "structure.out_ancestry_proportion_raceDemo_compare"),
#                            header = TRUE) |>
#     rename(BrNum = id)


dim(ancestry_tab)
# [1] 3134    4

metadata_visium_list <- read_csv(file = here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_visium_list.csv")) |> 
    left_join(ancestry_tab) |>
    mutate(BrNum = fct_reorder(BrNum, Eur))

samples_ancestry <- metadata_visium_list |> select(BrNum, Afr, Eur)
save(samples_ancestry, file = here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.Rdata"))
write.csv(samples_ancestry, file = here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.csv"), row.names = FALSE)

## local
# load(here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.Rdata"))
metadata_visium_list <- read_csv(file = here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_visium_list.csv")) |>
    left_join(samples_ancestry) 

metadata_visium_list <- metadata_visium_list |>
    mutate(BrNum = fct_reorder(BrNum, Eur),
           APOE = gsub(", ", "/", APOE),
           Ancestry = gsub('CAUC', "EA", Ancestry))

metadata_visium_list |> dplyr::count(APOE, Ancestry)

#### plot ancestry fractions ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

ancestry_colors2 <- ancestry_colors
names(ancestry_colors2) <- c("Eur", "Afr")

sample_ancestry_long <- metadata_visium_list |> 
    select(BrNum, APOE, Ancestry, Afr, Eur) |>
    pivot_longer(!c(BrNum, APOE, Ancestry), names_to = "anc_type", values_to = "anc_prop")

any(is.na(sample_ancestry_long$anc_prop))

sample_APOE_bar <- sample_ancestry_long |>
    group_by(APOE) |>
    mutate(APOE_anno = paste(APOE, "n=", n()/2)) |>
    ggplot(aes(x = BrNum, y = anc_prop, fill = anc_type)) +
    geom_col() +
    facet_wrap(~APOE_anno, scales = "free_x", ncol = 1) +
    scale_fill_manual(values = ancestry_colors2) +
    theme_bw() 

ggsave(sample_APOE_bar, filename = here(plot_dir, "Ancestry_APOE_bar.png"), width = 8)

# sample_ancestry_long |>
#     ggplot(aes(x = BrNum, y = anc_prop, fill = anc_type)) +
#     geom_col() +
#     facet_grid(APOE~Ancestry)

sample_ancestry_bar <- sample_ancestry_long |> 
    group_by(APOE, Ancestry) |>
    mutate(APOE_anno = paste(APOE, Ancestry, "n=", n()/2)) |>
    ggplot(aes(x = BrNum, y = anc_prop, fill = anc_type)) +
    geom_col() +
    facet_wrap(~APOE_anno, scales = "free_x", ncol = 2)   +
    scale_fill_manual(values = ancestry_colors2) +
    theme_bw() 

ggsave(sample_ancestry_bar, filename = here(plot_dir, "Ancestry_cat_APOE_bar.png"), width = 8)

#### Age ####

metadata_visium_list |> 
    ggplot(aes(y = Age, x = APOE, fill = APOE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(color = Ancestry)) +
    scale_color_manual(values = ancestry_colors) +
    scale_fill_manual(values = APOE_genotype_colors) +
    theme_bw()
    

# slurmjobs::job_single('04_ancestry_check', create_shell = TRUE, memory = '5G', command = "Rscript 04_ancestry_check.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
