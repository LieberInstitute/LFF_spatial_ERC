
library("tidyverse")
library("readxl")
library("here")

metadata_sn <- read_csv(here("processed-data", "00_project_prep", "02_get_online_metadata","metadata_sn_plan.csv"))
metadata_sn

dim(metadata_sn)
# [1] 31 16

## information about seqencing 
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
# [1] 31 26

write_csv(seq_info, here("processed-data", "04_snRNA-seq", "erc_sn_seq_info.csv"))

summary(seq_info)
# chromium_id           Tissue             BrNum           Experimental Round Sequencing Round   PI/NeuN         
# Length:14          Length:14          Length:14          Min.   :1.000      Min.   :1        Length:14         
# Class :character   Class :character   Class :character   1st Qu.:2.000      1st Qu.:1        Class :character  
# Mode  :character   Mode  :character   Mode  :character   Median :3.000      Median :1        Mode  :character  
# Mean   :2.643      Mean   :1                          
# 3rd Qu.:3.750      3rd Qu.:1                          
# Max.   :4.000      Max.   :1                          
# NA's   :2                          
#  Nuclei Sorted   Nuclei Targeted cDNA Amp Cycle Agilent [cDNA] pg/ul cDNA Dilution Factor Final [cDNA] pg/ul
#  Min.   :14000   Min.   :9000    Min.   :11     Min.   : 435.4       Min.   :1.000        Min.   : 435.4    
#  1st Qu.:14000   1st Qu.:9000    1st Qu.:11     1st Qu.:1158.1       1st Qu.:1.000        1st Qu.:1412.5    
#  Median :14000   Median :9000    Median :11     Median :1472.8       Median :1.000        Median :1812.1    
#  Mean   :14000   Mean   :9000    Mean   :11     Mean   :1635.9       Mean   :1.214        Mean   :1863.3    
#  3rd Qu.:14000   3rd Qu.:9000    3rd Qu.:11     3rd Qu.:2003.8       3rd Qu.:1.000        3rd Qu.:2263.9    
#  Max.   :14000   Max.   :9000    Max.   :11     Max.   :3304.7       Max.   :2.000        Max.   :3304.7    
#                                                                                                             
#  Total cDNA ng      cDNA Input       SI cycles     Ave frag length Agilent [lib] pg/ul Lib Dilution Factor
#  Min.   : 17.42   Min.   : 4.354   Min.   :14.00   Min.   :438.0   Min.   :1091        Min.   : 4.000     
#  1st Qu.: 56.50   1st Qu.:14.125   1st Qu.:14.00   1st Qu.:453.2   1st Qu.:2000        1st Qu.: 6.000     
#  Median : 72.48   Median :18.121   Median :14.00   Median :460.5   Median :2143        Median : 8.500     
#  Mean   : 74.53   Mean   :18.633   Mean   :14.57   Mean   :470.4   Mean   :2096        Mean   : 9.964     
#  3rd Qu.: 90.56   3rd Qu.:22.639   3rd Qu.:15.00   3rd Qu.:486.5   3rd Qu.:2431        3rd Qu.:12.250     
#  Max.   :132.19   Max.   :33.047   Max.   :16.00   Max.   :523.0   Max.   :2569        Max.   :24.000     
#                                                                                                           
#  Final [lib] pg/ul  index_name         index(i7)         index2_workflow_a(i5) index2_workflow_b(i5)
#  Min.   : 2355     Length:14          Length:14          Length:14             Length:14            
#  1st Qu.: 9228     Class :character   Class :character   Class :character      Class :character     
#  Median :14837     Mode  :character   Mode  :character   Mode  :character      Mode  :character     
#  Mean   :17912                                                                                      
#  3rd Qu.:21473                                                                                      
#  Max.   :47348                                                                                      
#                                                                                                     
#  Est Read Pairs (million)    Notes              round            sample_id        
#  Min.   :420              Length:14          Length:14          Length:14         
#  1st Qu.:420              Class :character   Class :character   Class :character  
#  Median :420              Mode  :character   Mode  :character   Mode  :character  
#  Mean   :420                                                                      
#  3rd Qu.:420                                                                      
#  Max.   :420                              
# # 

table(seq_info$round)

(cell_ranger_out_paths <- list.files(here("processed-data", "03_cellranger")))
# # [1] "10c_ERC_SVB" "11c_ERC_SVB" "12c_ERC_SVB" "13c_ERC_SVB" "14c_ERC_SVB" "15c_ERC_SVB" "16c_ERC_SVB" "17c_ERC_SVB"
# # [9] "18c_ERC_SVB" "19c_ERC_SVB" "1c_ERC_SVB"  "20c_ERC_SVB" "21c_ERC_SVB" "22c_ERC_SVB" "23c_ERC_SVB" "24c_ERC_SVB"
# # [17] "25c_ERC_SVB" "26c_ERC_SVB" "27c_ERC_SVB" "28c_ERC_SVB" "29c_ERC_SVB" "2c_ERC_SVB"  "30c_ERC_SVB" "31c_ERC_SVB"
# # [25] "3c_ERC_SVB"  "4c_ERC_SVB"  "5c_ERC_SVB"  "6c_ERC_SVB"  "7c_ERC_SVB"  "8c_ERC_SVB"  "9c_ERC_SVB"
# 
# length(cell_ranger_out_paths)
# # [1] 31

## Samples are all present in chromium files, and metadata
setequal(cell_ranger_out_paths, seq_info$chromium_id)
setdiff(cell_ranger_out_paths, seq_info$chromium_id)
#15-31
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

