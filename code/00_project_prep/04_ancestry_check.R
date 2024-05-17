
library(tidyverse)
library(here)
library(sessioninfo)

plot_dir <- here("plots", "00_project_prep", "04_ancestry_check")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

ancestry_tab <- read.table(here("raw-data", "ancestry", "structure.out_ancestry_proportion_raceDemo_compare"),
                           header = TRUE) |>
    rename(BrNum = id)
head(ancestry_tab)
# id   Afr   Eur Race
# 1 Br1857 0.000 1.000 CAUC
# 2 Br1306 0.791 0.209   AA
# 3 Br6171 0.000 1.000 CAUC
# 4 Br8404 0.000 1.000 CAUC
# 5 Br1802 0.854 0.146   AA
# 6 Br6345 0.735 0.265   AA

dim(ancestry_tab)
# [1] 2743    4

metadata_visium_list <- read_csv(file = here("processed-data", "00_project_prep", "02_get_online_metadata", "metadata_visium_list.csv")) |> 
    left_join(ancestry_tab) |>
    mutate(BrNum = fct_reorder(BrNum, Eur))

samples_ancestry <- metadata_visium_list |> select(BrNum, Afr, Eur)
save(samples_ancestry, file = here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.Rdata"))
write.csv(samples_ancestry, file = here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.csv"), row.names = FALSE)

sample_ancestry_long <- metadata_visium_list |> 
    select(BrNum, APOE, Ancestry, Afr, Eur) |>
    pivot_longer(!c(BrNum, APOE, Ancestry), names_to = "anc_type", values_to = "anc_prop")

any(is.na(sample_ancestry_long$anc_prop))

sample_APOE_bar <- sample_ancestry_long |>
    group_by(APOE) |>
    mutate(APOE_anno = paste(APOE, "n=", n()/2)) |>
    ggplot(aes(x = BrNum, y = anc_prop, fill = anc_type)) +
    geom_col() +
    facet_wrap(~APOE_anno, scales = "free_x", ncol = 1)

ggsave(sample_APOE_bar, filename = here(plot_dir, "Ancestry_APOE_bar.png"))

sample_ancestry_long
    ggplot(aes(x = BrNum, y = anc_prop, fill = anc_type)) +
    geom_col() +
    facet_grid(APOE~Ancestry)

sample_ancestry_bar <- sample_ancestry_long |> 
    group_by(APOE, Ancestry) |>
    mutate(APOE_anno = paste(APOE, Ancestry, "n=", n()/2)) |>
    ggplot(aes(x = BrNum, y = anc_prop, fill = anc_type)) +
    geom_col() +
    facet_wrap(~APOE_anno, scales = "free_x", ncol = 2)

ggsave(sample_ancestry_bar, filename = here(plot_dir, "Ancestry_cat_APOE_bar.png"))


