## Louise Huuki-Myers, April 2026
## Check polygenic risk score for dononrs

#### set up ####
library("tidyverse")
library("here")

data_dir <- here("processed-data", "00_project_prep", "12_polygenic_risk_score")
if(!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("plots", "00_project_prep", "12_polygenic_risk_score")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

#### load data ####

donor_flag <- read.csv(here("processed-data", "00_project_prep", "ERC_donor_notes.csv"))

donor_info <- read.csv(here("processed-data", "00_project_prep", "05_pathology", "sample_taupathy.csv")) |>
    filter(BrNum != "Br1289") |>
    mutate(donor_flag = BrNum %in% donor_flag$BrNum,
           APOE_carrier_Anc = paste(APOE_carrier, Ancestry))

PRS_scores <- read_csv(here(data_dir, "PRS_scores_by_threshold_APOEexcluded.csv")) |>
    pivot_longer(!IID, names_to = "cutoff", values_to = "PRS") |>
    rename(BrNum = IID) |>
    mutate(include_APOE = FALSE) |>
    bind_rows(read_csv(here(data_dir, "PRS_scores_by_threshold.csv")) |>
                  pivot_longer(!IID, names_to = "cutoff", values_to = "PRS") |>
                  rename(BrNum = IID) |>
                  mutate(include_APOE = TRUE)) |>
    mutate(cutoff_n = as.numeric(gsub("p cutoff=", "", cutoff)),
           cutoff = fct_reorder(gsub(" cutoff", "", cutoff), cutoff_n),
           include_APOE_anno = fct_rev(ifelse(include_APOE, "Include APOE", "Exclude APOE"))) |>
    right_join(donor_info, by = "BrNum")

PRS_scores |> count(cutoff)

PRS_scores_summary <- PRS_scores |> 
    group_by(cutoff, APOE_carrier, Ancestry, include_APOE) |> 
    summarize(min_PRS = min(PRS), median_PRS = median(PRS), max_PRS = max(PRS))

write_csv(PRS_scores_summary, file = here(data_dir, "ERC_PRS_scores_summary.csv"))

#### plots ####

prs_carrier_boxplot <- PRS_scores |> 
    ggplot(aes(x = cutoff, y = PRS, color = APOE_carrier)) +
    geom_boxplot() +
    # geom_jitter(width = .2) +
    facet_wrap(~include_APOE_anno, ncol = 1) +
    scale_color_manual(values = APOE_carrier_colors) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(prs_carrier_boxplot, filename = here(plot_dir, "PRS_carrier_boxplot.png"), width = 8, height = 6)

prs_carrier_boxplot_free <- PRS_scores |> 
    ggplot(aes(x = cutoff, y = PRS, color = APOE_carrier)) +
    geom_boxplot() +
    # geom_jitter(width = .2) +
    facet_wrap(~include_APOE_anno, scales = "free_y", ncol = 1) +
    scale_color_manual(values = APOE_carrier_colors) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw()  +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(prs_carrier_boxplot_free, filename = here(plot_dir, "PRS_carrier_boxplot_free.png"), width = 8, height = 6)

prs_carrier_boxplot_select <- PRS_scores |> 
    filter(cutoff %in% c("p=1e-08", "p=1e-05", "p=0.001", "p=0.05", "p=1")) |>
    ggplot(aes(x = cutoff, y = PRS, color = APOE_carrier)) +
    geom_boxplot() +
    # geom_jitter(position = ggplot2::position_jitter(seed = 1, width = .2)) +
    # ggrepel::geom_text_repel(aes(label = BrNum), position = ggplot2::position_jitter(seed = 1, width = .2), size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(include_APOE_anno~Ancestry) +
    scale_color_manual(values = APOE_carrier_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(prs_carrier_boxplot_select, filename = here(plot_dir, "PRS_carrier_boxplot_select.png"), width = 6, height = 4)

prs_carrier_boxplot05 <- PRS_scores |> 
    filter(cutoff == "p=0.05") |>
    ggplot(aes(x = APOE_carrier, y = PRS, color = APOE_carrier)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = ggplot2::position_jitter(seed = 1, width = .2)) +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = BrNum), position = ggplot2::position_jitter(seed = 1, width = .2), size = 2) +
    facet_grid(cutoff~include_APOE_anno) +
    scale_color_manual(values = APOE_carrier_colors) +
    labs(x = "APOE_carrier") +
    theme_bw() +
    theme(legend.position = "None")

ggsave(prs_carrier_boxplot05, filename = here(plot_dir, "PRS_carrier_boxplot_p05.png"), width = 4, height = 4)

#### Ancestry ####

prs_carrier_boxplot_anc <- PRS_scores |> 
    filter(cutoff == c("p=1e-08", "p=0.05","p=1")) |>
    # filter(cutoff %in% c("p=1e-08", "p=1e-05", "p=0.001", "p=0.05", "p=1")) |>
    ggplot(aes(y = PRS, x = APOE_carrier_Anc, fill = Ancestry)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = APOE_carrier), width = .2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(cutoff~include_APOE_anno, scales = "free_y") +
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = APOE_carrier_colors) +
    labs(x = "APOE_carrier + Anc") +
    theme_bw()
    

ggsave(prs_carrier_boxplot_anc, filename = here(plot_dir, "PRS_carrier_boxplot_Anc.png"))


prs_carrier_boxplot_flag <- PRS_scores |> 
    filter(cutoff == c("p=1e-08", "p=0.05","p=1")) |>
    # filter(cutoff %in% c("p=1e-08", "p=1e-05", "p=0.001", "p=0.05", "p=1")) |>
    ggplot(aes(y = PRS, x = APOE_carrier_Anc, fill = Ancestry)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = donor_flag), position = ggplot2::position_jitter(seed = 1, width = .2)) +
    ggrepel::geom_text_repel(aes(label = BrNum, color = donor_flag), position = ggplot2::position_jitter(seed = 1, width = .2), size = 2) +
    facet_grid(cutoff~include_APOE_anno, scales = "free_y") +
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = list(`TRUE` = "red", `FALSE` ="black")) +
    labs(x = "APOE_carrier + Anc") +
    theme_bw()
    

ggsave(prs_carrier_boxplot_flag, filename = here(plot_dir, "PRS_carrier_boxplot_flag.png"), height =12)

