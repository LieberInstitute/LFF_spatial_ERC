
library("tidyverse")
library("readxl")
library("here")

metadata_sn <- read.csv(here("processed-data", "00_project_prep", "02_get_online_metadata","metadata_sn_plan.csv"))
metadata_sn

dim(metadata_sn)
# [1] 31 16

## information about seqencing 
seq_info_fn <- list.files(here("raw-data", "sample_info"), full.names = TRUE)
names(seq_info_fn) <- c("round2", "round1")

seq_info <- map2(seq_info_fn, names(seq_info_fn), ~read_excel(.x, sheet = "Summary", range = "Summary!A1:W17") |>
                        mutate(Ct = as.double(Ct),
                               round = .y))

seq_info <- do.call("bind_rows", seq_info) |>
    filter(Tissue == "ERC") |>
    mutate(BrNum = gsub("\\(re-dis\\)","", gsub("Hs_", "", Brain)),
           sample_id = BrNum,
           chromium_id = gsub("_HRD", "", `Sample #`))


dim(seq_info)
# [1] 31 26

write_csv(seq_info, here("processed-data", "04_snRNA-seq", "erc_sn_seq_info.csv"))

summary(seq_info)
# Sample #            Tissue             Brain             Slide #            Array #                Ct       
# Length:31          Length:31          Length:31          Length:31          Length:31          Min.   :15.35  
# Class :character   Class :character   Class :character   Class :character   Class :character   1st Qu.:15.91  
# Mode  :character   Mode  :character   Mode  :character   Mode  :character   Mode  :character   Median :16.30  
# Mean   :16.39  
# 3rd Qu.:16.73  
# Max.   :18.14  
# NA's   :4      
#  cDNA Amp Cycle  Agilent [cDNA] pg/ul Dilution Factor...9 Final [cDNA] pg/ul Total cDNA ng yield   cDNA Input    
#  Min.   :16.00   Min.   : 874.5       Min.   :1.000       Min.   : 4722      Min.   :188.9       Min.   : 47.22  
#  1st Qu.:16.00   1st Qu.:1887.1       1st Qu.:3.000       1st Qu.: 8020      1st Qu.:320.8       1st Qu.: 80.20  
#  Median :17.00   Median :3288.7       Median :3.000       Median : 9875      Median :395.0       Median : 98.75  
#  Mean   :16.55   Mean   :3247.1       Mean   :4.161       Mean   :11051      Mean   :442.0       Mean   :110.51  
#  3rd Qu.:17.00   3rd Qu.:3947.0       3rd Qu.:5.000       3rd Qu.:13469      3rd Qu.:538.8       3rd Qu.:134.69  
#  Max.   :18.00   Max.   :7794.9       Max.   :9.000       Max.   :23385      Max.   :935.4       Max.   :233.85  
#                                                                                                                  
#    SI cycles     Ave frag length Agilent [lib] pg/ul Dilution Factor...16 Final [lib] pg/ul  index_name       
#  Min.   :15.00   Min.   :415.0   Min.   :1035        Min.   : 1.000       Min.   : 6208     Length:31         
#  1st Qu.:15.00   1st Qu.:443.0   1st Qu.:1326        1st Qu.: 4.000       1st Qu.:10998     Class :character  
#  Median :16.00   Median :453.0   Median :2046        Median : 7.000       Median :15615     Mode  :character  
#  Mean   :16.06   Mean   :453.3   Mean   :2151        Mean   : 7.032       Mean   :17833                       
#  3rd Qu.:17.00   3rd Qu.:465.5   3rd Qu.:2766        3rd Qu.: 9.500       3rd Qu.:22751                       
#  Max.   :17.00   Max.   :501.0   Max.   :3568        Max.   :23.000       Max.   :50161                       
#                  NA's   :8       NA's   :8                                NA's   :8                           
# index(i7)         index2_workflow_a(i5) index2_workflow_b(i5) % Coverage Array Est Read Pairs         round          
# Length:31          Length:31             Length:31             Min.   :45.00    Min.   :135000000   Length:31         
# Class :character   Class :character      Class :character      1st Qu.:67.50    1st Qu.:202500000   Class :character  
# Mode  :character   Mode  :character      Mode  :character      Median :75.00    Median :225000000   Mode  :character  
# Mean   :75.32    Mean   :225967742                     
# 3rd Qu.:85.00    3rd Qu.:255000000                     
# Max.   :98.00    Max.   :294000000                     
# 

# cn <- map(sample_info, colnames)
# table(unlist(cn))

# table(seq_info$round)
# round1 round2 
# 16     15

(cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger")))
# [1] "10c_ERC_SVB" "11c_ERC_SVB" "12c_ERC_SVB" "13c_ERC_SVB" "14c_ERC_SVB" "15c_ERC_SVB" "16c_ERC_SVB" "17c_ERC_SVB"
# [9] "18c_ERC_SVB" "19c_ERC_SVB" "1c_ERC_SVB"  "20c_ERC_SVB" "21c_ERC_SVB" "22c_ERC_SVB" "23c_ERC_SVB" "24c_ERC_SVB"
# [17] "25c_ERC_SVB" "26c_ERC_SVB" "27c_ERC_SVB" "28c_ERC_SVB" "29c_ERC_SVB" "2c_ERC_SVB"  "30c_ERC_SVB" "31c_ERC_SVB"
# [25] "3c_ERC_SVB"  "4c_ERC_SVB"  "5c_ERC_SVB"  "6c_ERC_SVB"  "7c_ERC_SVB"  "8c_ERC_SVB"  "9c_ERC_SVB"

length(cell_ranger_out_paths)
# [1] 31

cell_ranger_chrID <- gsub("c", "v", gsub("_SVB", "", cell_ranger_out_paths))

## Samples are all present in chromium files, and metadata
setequal(cell_ranger_chrID, seq_info$chromium_id)
setequal(metadata_sn$BrNum, seq_info$BrNum)


sample_info <- metadata_sn |>
    select(BrNum, Genotype, Age, Sex, Ancestry, Diagnosis, Rin, APOE) |>
    left_join(seq_info |>
                  select(sample_id, BrNum, chromium_id, round, Tissue, 
                         index_name, Array = `Array #`, Slide = `Slide #`)) 

sample_info$cell_ranger_path  <- cell_ranger_out_paths[order(as.integer(gsub("c_ERC_SVB","", cell_ranger_out_paths)))]

write_csv(sample_info, here("processed-data", "04_snRNA-seq", "erc_sn_sample_info.csv"))


all(file.exists(here("processed-data", "03_cellranger", sample_info$cell_ranger_path, "outs", "raw_feature_bc_matrix")))

## Create slurm job for droplet scores setup
slurmjobs::job_loop(
    loops = list(sample_id = cell_ranger_out_paths),
    name = "01_get_droplet_scores",
    create_shell = TRUE
)

