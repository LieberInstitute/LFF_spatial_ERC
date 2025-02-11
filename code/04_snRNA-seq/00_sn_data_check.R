# Aug, 2024 - Louise Huuki-Myers
#  Prep snRNA-seq data, Build metadata and collect sample info

library("tidyverse")
library("readxl")
library("here")

#### Build Metadata ####
## load ancestry data
load(here("processed-data","00_project_prep", "04_ancestry_check", "sample_ancestry.Rdata"), verbose = TRUE)
# samples_ancestry

# Info from Google doc, donor phenotype data
metadata_sn <- read_csv(here("processed-data", "00_project_prep", "02_get_online_metadata","metadata_sn_plan.csv")) |>
    mutate(APOE = gsub(", ", "/", APOE),
           Ancestry = gsub('CAUC', "EA", Ancestry),
           APOE_carrier = ifelse(grepl("E2", APOE), "E2+", "E4+")) |>
    left_join(samples_ancestry) |>
    rename(Anc_Afr = Afr, Anc_Eur = Eur)

metadata_sn 
metadata_sn  |> count(APOE, APOE_carrier) 


dim(metadata_sn)
# [1] 31 16

# Info from SVB, all 31 BrNum/chomiumIDs
erc_chrom_summary <- read_excel(here("processed-data", "04_snRNA-seq", "ERC_Chromium_summary_all.xlsx"), n_max = 39) |>
    filter(!is.na(`Brain #`)) |>
    rename(chromium_id = `Sample name`,
           BrNum = `Brain #`) |>
    mutate(BrNum = gsub(" ", "", BrNum))

erc_chrom_summary |> count(Round)

## information about sequencing - but only first 14 samples
seq_info_fn <- list.files(here("raw-data", "sample_info_snRNAseq"), full.names = TRUE)
names(seq_info_fn) <- gsub("^.*?ERC_chromium_(round_[1-4])_sequencing_summary.xlsx", "\\1", basename(seq_info_fn))

seq_info <- map2(seq_info_fn, names(seq_info_fn), ~read_excel(.x, sheet = "Summary") |>
                        mutate(round = .y))

seq_info <- do.call("bind_rows", seq_info) |>
    rename(chromium_id = `Sample #`,
           BrNum = Brain) |>
    mutate(sample_id = BrNum)

seq_info |> count(`Experimental Round`, round)

dim(seq_info)
# [1] 14 27

# write_csv(seq_info, here("processed-data", "04_snRNA-seq", "erc_sn_seq_info.csv"))

## combine info from erc_chrom_summary + seq_info
seq_info2 <- erc_chrom_summary |>
    select(BrNum, chromium_id, exp_round = Round) |> 
    left_join(seq_info |> select(BrNum, chromium_id, seq_round = `Sequencing Round`, exp_round = `Experimental Round`)) |>
    replace_na(list(seq_round = 2))

seq_info2 |> count(seq_round)
# seq_round     n
# <dbl> <int>
# 1         1    12
# 2         2    19
seq_info2 |> count(exp_round)
# exp_round     n
# <dbl> <int>
# 1         1     3
# 2         2     3
# 3         3     4
# 4         4     4
# 5         5     4
# 6         6     4
# 7         7     4
# 8         8     5
seq_info2 |> count(seq_round, exp_round)

(cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger")))
# # [1] "10c_ERC_SVB" "11c_ERC_SVB" "12c_ERC_SVB" "13c_ERC_SVB" "14c_ERC_SVB" "15c_ERC_SVB" "16c_ERC_SVB" "17c_ERC_SVB"
# # [9] "18c_ERC_SVB" "19c_ERC_SVB" "1c_ERC_SVB"  "20c_ERC_SVB" "21c_ERC_SVB" "22c_ERC_SVB" "23c_ERC_SVB" "24c_ERC_SVB"
# # [17] "25c_ERC_SVB" "26c_ERC_SVB" "27c_ERC_SVB" "28c_ERC_SVB" "29c_ERC_SVB" "2c_ERC_SVB"  "30c_ERC_SVB" "31c_ERC_SVB"
# # [25] "3c_ERC_SVB"  "4c_ERC_SVB"  "5c_ERC_SVB"  "6c_ERC_SVB"  "7c_ERC_SVB"  "8c_ERC_SVB"  "9c_ERC_SVB"
# 
# length(cell_ranger_out_paths)
# # [1] 31

## Samples are all present in chromium files, and metadata
setequal(cell_ranger_out_paths, seq_info2$chromium_id)
setequal(metadata_sn$BrNum, seq_info2$BrNum)

sample_info <- metadata_sn |>
    select(BrNum, APOE, APOE_carrier, Ancestry, Sex, Age, Diagnosis, Rin, Anc_Afr, Anc_Eur) |>
    left_join(seq_info2) |>
    mutate(chromium_index = as.integer(gsub("c_ERC_SVB","",chromium_id)),
           sample_id = BrNum,
           exp_round = paste0("round", exp_round),
           seq_round = paste0("round", seq_round)
           ) |>
    arrange(chromium_index) |>
    select(sample_id, chromium_id, BrNum, everything())

write_csv(sample_info, here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))

## Create slurm job for droplet scores setup
slurmjobs::job_loop(
    loops = list(sample_id = cell_ranger_out_paths),
    name = "01_get_droplet_scores",
    create_shell = TRUE
)

#### manual knee values ####
manual_knee <- c("8c_ERC_SVB", "5c_ERC_SVB", "24c_ERC_SVB", "27c_ERC_SVB", "28c_ERC_SVB", "29c_ERC_SVB")
manual_knee_lgl <- cell_ranger_out_paths %in% manual_knee
names(manual_knee_lgl) <- cell_ranger_out_paths

paste(manual_knee_lgl, collapse = " ")


slurmjobs::array_submit(
    name = "01_get_droplet_scores",
    task_ids = c(which(manual_knee_lgl) -1),
    submit = TRUE
)

