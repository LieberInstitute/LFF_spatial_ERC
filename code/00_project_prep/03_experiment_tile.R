
library("tidyverse")
library("googlesheets4")
library("sessioninfo")
library("here")

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

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

metadata_plan |> dplyr::count(assay, complete)
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
    
ggsave(experiment_tile, filename = here(plot_dir, "experiment_tile.png"), height = 4, width = 7)

#### Donor Info ####

sample_info <- read_csv(here("processed-data", "02_build_spe", "sample_info.csv"))
sample_info$APOE_carrier = ifelse(grepl("E2", sample_info$APOE), "E2+", "E4+")

table(sample_info$Sex)
# F  M 
# 9 22

table(sample_info$Ancestry)
# AA EA 
# 17 14 

table(sample_info$APOE_carrier)
# E2+ E4+ 
# 14  17

table(sample_info$APOE)
# E2/E2 E2/E3 E3/E4 E4/E4 
# 6     8    10     7 

summary(sample_info$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 29.95   47.61   51.63   52.36   59.92   68.38

# genotype tile plot
APOE_ancestry_tile <- sample_info |> 
    dplyr::count(Ancestry, APOE) |>
    ggplot(aes(x = APOE, y = Ancestry, 
               fill = APOE)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), color = "black") +
    scale_fill_manual(values = APOE_genotype_colors) +
    # scale_color_manual(values = ancestry_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(APOE_ancestry_tile, filename = here(plot_dir, "LFF_ERC_APOE_ancestry_tileplot.png"), height = 2, width = 3)

## carrier tile plot
APOE_carrier_count <- sample_info |> 
    group_by(Ancestry, APOE_carrier) |>
    summarise(n = n(), 
              n_M = sum(Sex == "M"),
              n_F = sum(Sex == "F")) |>
    mutate(anno = sprintf("%i\n(m:%i, f:%i)", n, n_M, n_F))

APOE_carrier_ancestry_tile <- APOE_carrier_count |>
    ggplot(aes(x = APOE_carrier, y = Ancestry, 
               fill = APOE_carrier)) +
    geom_tile(color = "white") +
    geom_text(aes(label = anno), color = "black") +
    scale_fill_manual(values = APOE_carrier_colors) +
    # scale_color_manual(values = ancestry_colors) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(APOE_carrier_ancestry_tile, filename = here(plot_dir, "LFF_ERC_APOE_carrier_ancestry_tileplot.png"), height = 2, width = 3)

## tile plot split by sex
APOE_ancestry_sex_tile <- sample_info |> 
    dplyr::count(Ancestry, APOE, Sex) |>
    ggplot(aes(x = APOE, y = Ancestry, 
               fill = APOE)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), color = "black") +
    scale_fill_manual(values = APOE_genotype_colors) +
    facet_wrap(~Sex) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(APOE_ancestry_sex_tile, filename = here(plot_dir, "LFF_ERC_APOE_ancestry_sex_tileplot.png"), height = 2, width = 6)

## sex barplot
sex_barplot <- sample_info |>
    ggplot(aes(x = Sex, fill = Sex)) +
    geom_bar() +
    geom_text(stat='count',aes(label = after_stat(count)), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = sex_colors) +
    theme_bw() +
    theme(legend.position= "None")

ggsave(sex_barplot, filename = here(plot_dir, "sex_barplot.png"), width = 2, height = 2)


## Age boxplot + genotype + sex

age_boxplot <- sample_info |>
    ggplot(aes(x = APOE, y = Age, color = Sex)) +
    geom_boxplot() +
    facet_wrap(~Ancestry) +
    theme_bw()

ggsave(age_boxplot, filename = here(plot_dir, "age_boxplot.png"))


# slurmjobs::job_single('03_experiment_tile', create_shell = TRUE, memory = '5G', command = "Rscript 03_experiment_tile.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

